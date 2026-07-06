"""
run_visco.py -- viscoelastic curing reference model, CalculiX pipeline.

Steps:
  1. (optional) regenerate the CalculiX deck:   python convert_visco_inp.py
  2. compile libviscmod.dll from
     Abaqus-Viscoelastic-Curing-Subroutine_calculix.for + umat_shim.f
  3. run ccx_ext.exe -i curing_visco_ccx
  4. parse the .dat element histories (16-element through-thickness
     column at the plate centre) and plot transverse in-ply stress vs
     time against the digitized Abaqus reference curve.

Usage:
  python run_visco.py               # solve + plot (compile only if DLL missing)
  python run_visco.py --compile     # force recompiling libviscmod.dll
  python run_visco.py --plot-only   # skip DLL build + solve
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
CCX = Path(r"C:\YVAN\CODE\CalculiX\src\ccx_ext.exe")
GFORTRAN = Path(r"C:\msys64\ucrt64\bin\gfortran.exe")
JOB = "curing_visco_ccx"
SRC_FOR = "Abaqus-Viscoelastic-Curing-Subroutine_calculix.for"

# through-thickness output column (bottom -> top), see convert_visco_inp.py
OUTCOL = [7500, 7475, 7450, 7425, 2401, 2426, 2451, 2476,
          10000, 9975, 9950, 9925, 2525, 2550, 2575, 2600]
# ply of each column element ('0' = fiber x, '90' = fiber y)
PLY = ["0"] * 4 + ["90"] * 8 + ["0"] * 4

# digitized "Reference solution" curve of the Chile benchmark plot
# (stress [MPa] vs time [min], plateau ~17.5, final ~41.5)
REF_T = np.array([0, 40, 60, 90, 97, 105, 113, 120, 125, 132, 140,
                  150, 160, 175, 200, 240, 260, 280, 300])
REF_S = np.array([0, 0.4, 0.6, 1.0, 1.2, 0.9, 0.6, 1.0, 2.2, 5.5, 9.5,
                  12.8, 14.5, 15.8, 17.0, 17.5, 25.5, 33.5, 41.5])


def build_dll():
    cmd = [str(GFORTRAN), "-shared", "-O2", "-ffixed-form",
           "-ffixed-line-length-none", "-Wl,--export-all-symbols",
           "-J", "_visco_build", "-I", "_visco_build",
           SRC_FOR, "umat_shim.f", "-o", "libviscmod.dll"]
    (HERE / "_visco_build").mkdir(exist_ok=True)
    print("[build]", " ".join(cmd))
    subprocess.run(cmd, cwd=HERE, check=True)


def run_ccx():
    print(f"[ccx] {CCX} -i {JOB}")
    r = subprocess.run([str(CCX), "-i", JOB], cwd=HERE)
    print("[ccx] exit", r.returncode)
    if r.returncode:
        sys.exit(r.returncode)


def parse_dat():
    """-> times [s], stress{eid: (n_t, 6)}, alpha{eid: (n_t,)},
    temp{eid: (n_t,)} [K]"""
    stress, sdv = {}, {}
    hdr = re.compile(
        r"(stresses|internal state variables).*time\s+([\d.E+\-]+)", re.I)
    mode, t = None, None
    with open(HERE / f"{JOB}.dat", encoding="latin-1") as f:
        for raw in f:
            line = raw.strip()
            m = hdr.match(line)
            if m:
                t = float(m.group(2))
                mode = "s" if m.group(1).startswith("stress") else "v"
                d = stress if mode == "s" else sdv
                d.setdefault(t, {})
                continue
            if mode and line and re.match(r"[a-zA-Z]", line):
                mode = None
                continue
            if mode and line:
                p = line.split()
                try:
                    eid = int(p[0])
                except ValueError:
                    continue
                # lines of oriented elements end with the orientation name
                d = stress if mode == "s" else sdv
                if mode == "s":
                    d[t].setdefault(eid, []).append(
                        [float(v) for v in p[2:8]])
                else:
                    # SDV1 = alpha, SDV10 = temperature [K]
                    d[t].setdefault(eid, []).append(
                        [float(p[2]), float(p[11])])
    times = sorted(stress)
    S = {e: np.array([np.mean(stress[t].get(e, [[np.nan] * 6]), axis=0)
                      for t in times]) for e in OUTCOL}
    sv = {e: np.array([np.mean(sdv[t].get(e, [[np.nan] * 2]), axis=0)
                       for t in times if t in sdv]) for e in OUTCOL}
    A = {e: sv[e][:, 0] for e in OUTCOL}
    TK = {e: sv[e][:, 1] for e in OUTCOL}
    return np.array(times), S, A, TK


def make_plot(times, S, A):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tmin = times / 60.0
    fig, (ax, axa) = plt.subplots(2, 1, figsize=(10, 8.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(REF_T, REF_S, "r--", lw=2.4,
            label="Chile benchmark reference solution (digitized)")

    # elements carry *ORIENTATION, so .dat stresses are in the local
    # material frame: transverse-to-fiber in-ply stress is local S22
    # for every ply. Show mid-thickness 90-ply and bottom/top 0-ply.
    picks = [(OUTCOL[6], "90", "mid 90-ply (z~10.3 mm)", "tab:blue"),
             (OUTCOL[1], "0", "bottom 0-ply (z~2.4 mm)", "tab:green"),
             (OUTCOL[14], "0", "top 0-ply (z~23 mm)", "tab:purple")]
    for eid, ply, label, col in picks:
        comp = 1   # local S22 = transverse in-ply
        ax.plot(tmin, S[eid][:, comp] / 1e6, color=col, lw=1.6,
                label=f"CalculiX {label}, transverse in-ply")
    for tt in (2184, 5784, 7248, 14448):
        ax.axvline(tt / 60, color="k", lw=0.5, alpha=0.35)
    ax.set_ylabel("transverse in-ply stress [MPa]")
    ax.axhline(0, color="k", lw=0.7)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_title("[0/90/0] laminate, two-dwell MRCC cycle -- viscoelastic DLL "
                 "(libviscmod.dll) via ccx_ext.exe", fontsize=10.5)

    amean = np.mean([A[e] for e in OUTCOL], axis=0)
    axa.plot(tmin[:len(amean)], amean, color="0.4", lw=1.5)
    axa.axhline(0.30, color="0.6", lw=0.8, ls=":")
    axa.text(2, 0.33, "gel", fontsize=8, color="0.4")
    axa.set_ylabel("alpha (column mean)")
    axa.set_xlabel("Time [min]")
    axa.grid(True, alpha=0.3)

    fig.suptitle("Viscoelastic curing reference model: CalculiX vs Abaqus",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    out = HERE / "curing_visco_vs_reference.png"
    fig.savefig(out, dpi=110)
    print("[plot] saved", out)


def main():
    if "--plot-only" not in sys.argv:
        dll = HERE / "libviscmod.dll"
        if "--compile" in sys.argv or not dll.exists():
            build_dll()
        else:
            print(f"[build] {dll.name} up to date? not checked -- reusing "
                  "existing DLL (pass --compile to force a rebuild)")
        run_ccx()
    times, S, A, _ = parse_dat()
    print(f"[dat] {len(times)} output frames, last t = {times[-1]:.0f} s")
    make_plot(times, S, A)


if __name__ == "__main__":
    main()
