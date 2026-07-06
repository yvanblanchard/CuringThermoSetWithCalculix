"""
compare_simplecure_visco.py -- run the 1-D analytical curing model
(simplecure_autoclave.py, SimpleCure port) with the SAME inputs as the
viscoelastic CalculiX FEA model (curing_visco_ccx.inp / libviscmod.dll)
and compare the degree-of-cure histories.

Same inputs means:
  - geometry      : 25.4 mm laminate, no tool (cycle applied on both faces,
                    h very large to mimic the FEA prescribed-temperature BC)
  - cure cycle    : two-dwell MRCC (25C -> 2.5 K/min -> 116C x 60 min ->
                    2.5 K/min -> 177C x 120 min); the final cool-down of the
                    FEA cycle is omitted (cure is complete before it, and the
                    1-D tool only supports heating ramps)
  - thermal props : rho_c = 1578 kg/m3, cp = 862 J/kgK, k_transverse =
                    0.4135 W/mK (deck values)
  - kinetics      : Lee-Loos 3501-6 exactly as in the UMAT
                    (Abaqus-Viscoelastic-Curing-Subroutine_calculix.for):
                       alpha <= 0.3 : dadt = (K1+K2*a)(1-a)(B-a)   [1/min]
                       alpha  > 0.3 : dadt = K3*(1-a)              [1/min]
                    injected in place of simplecure's dadt1 (neither of its
                    built-in models can represent Lee-Loos)
  - exotherm      : q = (1-Vf) * rho_r * Htot * dadt  (same form both codes)
  - Tg            : DiBenedetto, Tg0 = -20C, TgInf = 220C, lambda = 0.5

CalculiX side: alpha(t) = SDV1 histories of the through-thickness OUTCOL
column parsed from curing_visco_ccx.dat (requires a completed run).

Usage:  python compare_simplecure_visco.py
Output: simplecure_vs_visco.png + printed checkpoint table
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import simplecure_autoclave as sc
import run_visco as rv

# ---- UMAT constants (deck *USER MATERIAL, see convert_visco_inp.py) -------
A1, A2, A3 = 2.101e9, -2.014e9, 1.960e5      # [1/min]
DE1, DE2, DE3 = 8.07e4, 7.78e4, 5.66e4       # [J/mol]
B_CAT = 0.47
RHO_R, HTOT, VF = 1272., 4.736e5, 0.52
RHO_C, CP_C, K_TRANS = 1578., 862., 0.4135
TG0_C, TGINF_C, LAMBDA = -20., 220., 0.5


def lee_loos_dadt1(r, a1, e1, a2, e2, m, n1, n2, ad, ed, b, w, g,
                   tg, temp, conv):
    """Lee-Loos 3501-6 kinetics, dadt in 1/s (signature of sc.dadt1).

    Mapping: a1/e1 -> A1/DE1, a2/e2 -> A2/DE2, ad/ed -> A3/DE3, b -> B.
    """
    tk = temp + 273.15
    if 1.0 - conv <= 1e-3:
        return 0.0
    if conv <= 0.3:
        k1 = a1 * math.exp(-e1 / (r * tk))
        k2 = a2 * math.exp(-e2 / (r * tk))
        rate = (k1 + k2 * conv) * (1.0 - conv) * (b - conv)
    else:
        k3 = ad * math.exp(-ed / (r * tk))
        rate = k3 * (1.0 - conv)
    return max(rate, 0.0) / 60.0             # [1/min] -> [1/s]


def build_inputs():
    """Inputs grid with the FEA deck values (VBA-style Cells(row, col))."""
    rho_f = (RHO_C - (1 - VF) * RHO_R) / VF   # so composite density = deck's
    g = {
        (1, 2): 25.4,     # part thickness [mm]
        (2, 2): 0.0,      # tool thickness [mm] -- no tool in the FEA model
        (5, 2): 25.0,     # initial T [C]  (298.15 K)
        (6, 2): 1e-6,     # initial alpha (UMAT floor)
        (9, 2): 2.5, (10, 2): 116., (11, 2): 60.,    # ramp1, dwell1
        (12, 2): 2.5, (13, 2): 177., (14, 2): 120.,  # ramp2, dwell2
        (15, 2): 0.0,                                 # no 3rd ramp
        (20, 2): 1e5, (21, 2): 1e5,  # h top/bottom: mimic prescribed T
        (23, 2): 150.,               # Tg threshold for 'cure time' [C]
        # resin: density + constant cp/k giving the deck composite values
        (2, 5): RHO_R,
        (5, 5): CP_C,                # cpr = brcp = 862 (rom_cp=0 -> used as-is)
        (13, 5): K_TRANS,            # kr = dkr = 0.4135 (rom_k=0 -> as-is)
        (20, 5): VF,
        (2, 8): rho_f,               # fiber density (recovers rho_c = 1578)
        (13, 8): 0., (14, 8): 0., (15, 8): 0.,  # tool props (unused)
        (18, 8): 1,                  # kin_flag = 1 -> our injected Lee-Loos
        (19, 8): 0, (20, 8): 0,      # rom_cp / rom_k off (constants above)
        # kinetics model-1 cells -> Lee-Loos constants (see mapping above)
        (2, 11): A1, (3, 11): DE1, (4, 11): A2, (5, 11): DE2,
        (9, 11): A3, (10, 11): DE3, (11, 11): B_CAT,
        (14, 11): HTOT,
        (18, 11): TG0_C, (19, 11): TGINF_C, (20, 11): LAMBDA,
    }
    grid = {}
    for (row, col), v in g.items():
        grid.setdefault(row, {})[col] = v
    return sc.Inputs(grid)


def main():
    sc.dadt1 = lee_loos_dadt1                 # inject UMAT kinetics
    out = sc.run_simulation(build_inputs())
    sc.print_results(out["results"])

    t1d = np.array(out["time_min"])           # [min]
    x1d = np.array(out["x_positions"])        # [mm]
    a1d = np.array(out["degree_of_cure"])     # [time, node]
    T1d = np.array(out["temperature"])
    air = np.array(out["thermal_profile"])

    imid = int(np.argmin(np.abs(x1d - 12.7)))
    isurf = int(np.argmin(np.abs(x1d - 0.8)))  # ~ centroid of 1st FEA element

    # ---- CalculiX alpha + temperature histories ---------------------------
    times, S, A, TK = rv.parse_dat()
    tfe = times / 60.0
    acol = np.array([A[e] for e in rv.OUTCOL])          # [16, time]
    tcol = np.array([TK[e] for e in rv.OUTCOL]) - 273.15  # [C]
    zc = (np.arange(16) + 0.5) * 25.4 / 16.0            # element centroids
    afe_mid = acol[[7, 8]].mean(axis=0)                 # straddle z = 12.7
    afe_surf = acol[0]                                  # z ~ 0.79 mm
    Tfe_mid = tcol[[7, 8]].mean(axis=0)
    Tfe_surf = tcol[0]

    # ---- checkpoint table ---------------------------------------------------
    print("\n t[min]   alpha mid 1D / FEA      alpha surf 1D / FEA")
    for tm in (40, 60, 90, 110, 120, 130, 140, 160, 200, 240):
        i1 = np.argmin(np.abs(t1d - tm)); i2 = np.argmin(np.abs(tfe - tm))
        print(f" {tm:6.0f}   {a1d[i1, imid]:.3f} / {afe_mid[i2]:.3f}"
              f"            {a1d[i1, isurf]:.3f} / {afe_surf[i2]:.3f}")

    # ---- plot ---------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, axt) = plt.subplots(2, 1, figsize=(10, 8.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 2]})
    ax.plot(tfe, afe_mid, "b", lw=1.8, label="CalculiX FEA, mid-thickness")
    ax.plot(tfe, afe_surf, "b--", lw=1.4, label="CalculiX FEA, near surface")
    ax.plot(t1d, a1d[:, imid], "r", lw=1.8, label="SimpleCure 1D, mid-thickness")
    ax.plot(t1d, a1d[:, isurf], "r--", lw=1.4, label="SimpleCure 1D, near surface")
    ax.axhline(0.30, color="0.6", lw=0.8, ls=":")
    ax.text(2, 0.32, "gel / kinetics branch switch", fontsize=8, color="0.4")
    ax.set_ylabel("degree of cure alpha [-]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="lower right")
    ax.set_title("Lee-Loos 3501-6 kinetics: CalculiX viscoelastic FEA vs "
                 "SimpleCure 1D analytical model (same inputs)", fontsize=10.5)

    axt.plot(t1d, air, "k--", lw=1.2, label="autoclave air (cycle)")
    axt.plot(tfe, Tfe_mid, "b", lw=1.5, label="CalculiX FEA, mid-thickness")
    axt.plot(tfe, Tfe_surf, "b--", lw=1.1, label="CalculiX FEA, near surface")
    axt.plot(t1d, T1d[:, imid], "r", lw=1.5, label="SimpleCure 1D, mid-thickness")
    axt.set_ylabel("temperature [C]")
    axt.set_xlabel("Time [min]")
    axt.grid(True, alpha=0.3)
    axt.legend(fontsize=9, loc="lower right")

    fig.tight_layout()
    outpng = HERE / "simplecure_vs_visco.png"
    fig.savefig(outpng, dpi=110)
    print("\n[plot] saved", outpng)


if __name__ == "__main__":
    main()
