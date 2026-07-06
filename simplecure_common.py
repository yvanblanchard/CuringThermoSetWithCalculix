"""
SimpleCure -- shared physics and single-domain FE solver (Python port of the
VBA in SimpleCure_Infusion.xlsm / SimpleCure_Oven.xlsm / SimpleCure_RTM.xlsm).

These three workbooks model a single composite laminate (no separate tool) with a
1-D linear finite-element, theta-method (theta = 0.999) transient heat solve
coupled to resin cure kinetics and glass-transition temperature (Tg).  They
differ only in the boundary conditions applied to the two faces:

    Infusion : bottom face prescribed = thermal profile (heated mould);
               top face natural convection to ambient `too` (h from Inputs).
    Oven     : bottom face prescribed = thermal profile;
               top face convection to the oven air (sink = thermal profile).
    RTM      : both faces prescribed = thermal profile (closed heated mould),
               no convection.

The VBA helper subroutines dadt1 / dadt2 / kc / cpc are identical to the
autoclave workbook and are re-implemented here.  The custom inverse / mmult /
transpose routines are replaced by numpy.linalg.solve.

Faithful-porting note: the Infusion/Oven/RTM macros reference an *undeclared*
variable `dxp_last` inside the element-matrix assembly and the overshoot-location
code.  Under VBA's implicit typing that name is an empty Variant (= 0), so every
`If (dxp_last > 0)` branch there is dead code and the last partial element gets
no special treatment in the stiffness/mass matrices.  Only the separately-named
`dx_last` (used for nj(nn) and the node x-positions) is live.  This module
reproduces that behaviour exactly.
"""

import math
import os
import csv

import numpy as np

try:
    import openpyxl
except ImportError:
    openpyxl = None


# ---------------------------------------------------------------------------
# Cure kinetics / material property models (identical to the autoclave port)
# ---------------------------------------------------------------------------
def dadt1(r, a1, e1, a2, e2, m, n1, n2, ad, ed, b, w, g, tg, temp, conv):
    """Cure kinetics model 1 (diffusion-limited). Returns da/dt."""
    k1c = a1 * math.exp(-e1 / (r * (temp + 273.15)))
    k2c = a2 * math.exp(-e2 / (r * (temp + 273.15)))
    f = w * (temp - tg) + g
    if f == 0:
        k1, k2 = k1c, k2c
    else:
        kd = ad * math.exp(-ed / (r * (temp + 273.15))) * math.exp(-b / f)
        k1 = 0.0 if k1c == 0 else 1.0 / (1.0 / k1c + 1.0 / kd)
        k2 = 0.0 if k2c == 0 else 1.0 / (1.0 / k2c + 1.0 / kd)
    if (1 - conv) > 0.001:
        return k1 * (1 - conv) ** n1 + k2 * conv ** m * (1 - conv) ** n2
    return 0.0


def dadt2(r, a11, e11, a22, e22, m11, n11, m22, n22, d11, aco, act, temp, conv):
    """Cure kinetics model 2 (autocatalytic, two-stage). Returns da/dt."""
    k1 = a11 * math.exp(-e11 / (r * (temp + 273.15)))
    k2 = a22 * math.exp(-e22 / (r * (temp + 273.15)))
    if (1 - conv) > 0.001:
        return (k1 * conv ** m11 * (1 - conv) ** n11
                + (k2 * conv ** m22 * (1 - conv) ** n22)
                / (1 + math.exp(d11 * (conv - (aco + act * (temp + 273.15))))))
    return 0.0


def kc(rom_k, vf, akr, bkr, ckr, dkr, ekr, fkr, akf, bkf, temp, conv):
    """Composite thermal conductivity (Springer-Tsai rule of mixtures optional)."""
    kr = (akr * temp * conv + bkr * conv + ckr * temp + dkr
          + ekr * temp * conv ** 2 + fkr * conv ** 2)
    kf = akf * temp + bkf
    if rom_k > 0:
        return (vf * kr * (kf / kr - 1)
                + kr * (0.5 - kf / (2 * kr))
                + kr * (kf / kr - 1)
                * (vf ** 2 - vf + ((kf / kr + 1) ** 2) / ((2 * kf / kr - 2) ** 2)) ** 0.5)
    return kr


def cpc(rom_cp, tg, arcp, brcp, drcp, crcp, sigma, afcp, bfcp, wf, temp):
    """Composite specific heat (rule of mixtures optional)."""
    cpr = arcp * temp + brcp + drcp / (1 + math.exp(crcp * (temp - tg - sigma)))
    cpf = afcp * temp + bfcp
    if rom_cp > 0:
        return wf * cpf + (1 - wf) * cpr
    return cpr


# ---------------------------------------------------------------------------
# Inputs sheet accessor (VBA-style 1-based Cells(row, col))
# ---------------------------------------------------------------------------
class Inputs:
    def __init__(self, grid):
        self._grid = grid

    def cell(self, row, col):
        try:
            v = self._grid[row][col]
        except (IndexError, KeyError):
            v = None
        return 0.0 if v is None else v

    @classmethod
    def from_workbook(cls, path):
        if openpyxl is None:
            raise RuntimeError("openpyxl is required. Install with: pip install openpyxl")
        wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
        ws = wb["Inputs"]
        grid = {}
        for r in ws.iter_rows():
            for c in r:
                if c.value is not None:
                    grid.setdefault(c.row, {})[c.column] = c.value
        return cls(grid)


# ---------------------------------------------------------------------------
# Generic single-domain solver (Infusion / Oven / RTM)
# ---------------------------------------------------------------------------
def solve_single_domain(inp, *, dirichlet_top, convection_top, top_sink_mode,
                        threshold_row, series_name_profile):
    """Run a single-domain autoclave-family cure simulation.

    Parameters
    ----------
    inp : Inputs
    dirichlet_top : bool
        True  -> top face temperature prescribed = thermal profile (RTM).
        False -> top face is a free/convective boundary (Infusion, Oven).
    convection_top : bool
        True  -> apply surface convection at the top node (Infusion, Oven).
    top_sink_mode : {'ambient', 'profile', None}
        Convection sink temperature source: 'ambient' uses Too from Inputs
        (Infusion), 'profile' uses the current thermal profile (Oven), None
        means no convection (RTM).
    threshold_row : int
        Row in column B holding the fully-cured Tg threshold (23 for
        Infusion/Oven, 19 for RTM).
    series_name_profile : str
        Legend label for the driving-temperature series ("Mould Temperature",
        "Oven Temperature", ...).
    """
    Cells = inp.cell

    dx = 0.0005
    dt_nom = 30
    theta = 0.999

    # ---- single-domain mesh ----------------------------------------------
    ne_dx = int(math.floor(Cells(1, 2) / 1000 / dx))
    dx_last = (Cells(1, 2) / 1000) - (ne_dx * dx)
    ne = ne_dx if dx_last == 0 else ne_dx + 1
    nn = ne + 1

    # ---- resin properties -------------------------------------------------
    pr = Cells(2, 5)
    arcp, brcp, crcp, drcp, sigma = (Cells(4, 5), Cells(5, 5), Cells(6, 5),
                                     Cells(7, 5), Cells(8, 5))
    akr, bkr, ckr, dkr, ekr, fkr = (Cells(10, 5), Cells(11, 5), Cells(12, 5),
                                    Cells(13, 5), Cells(14, 5), Cells(15, 5))

    # ---- fibre properties -------------------------------------------------
    pf = Cells(2, 8)
    afcp, bfcp = Cells(4, 8), Cells(5, 8)
    akf, bkf = Cells(7, 8), Cells(8, 8)

    vf = Cells(20, 5)
    p = vf * pf + (1 - vf) * pr
    wf = vf * pf / p

    # ---- convection (Infusion/Oven only) ---------------------------------
    h = Cells(20, 2) if convection_top else 0.0
    too = Cells(21, 2)  # ambient sink (Infusion)

    # ---- initial conditions ----------------------------------------------
    tinit = Cells(5, 2)
    ainit = Cells(6, 2)

    # ---- thermal profile parameters --------------------------------------
    ramp1 = Cells(9, 2) / 60
    temp1 = Cells(10, 2)
    time1 = Cells(11, 2) * 60
    ramp2 = Cells(12, 2) / 60
    temp2 = Cells(13, 2)
    time2 = Cells(14, 2) * 60
    ramp3 = Cells(15, 2) / 60
    temp3 = Cells(16, 2)
    time3 = Cells(17, 2) * 60
    if ramp2 == 0:
        ramp3 = 0

    # ---- increment / time-step schedule ----------------------------------
    time_ramp1 = (temp1 - tinit) / ramp1
    inc_dt_ramp1 = int(math.floor(time_ramp1 / dt_nom))
    dt_last_ramp1 = time_ramp1 - inc_dt_ramp1 * dt_nom
    inc_ramp1 = inc_dt_ramp1 if dt_last_ramp1 == 0 else inc_dt_ramp1 + 1

    inc_dt_dwell1 = int(math.floor(time1 / dt_nom))
    dt_last_dwell1 = time1 - inc_dt_dwell1 * dt_nom
    inc_dwell1 = inc_dt_dwell1 if dt_last_dwell1 == 0 else inc_dt_dwell1 + 1

    if ramp2 > 0:
        time_ramp2 = (temp2 - temp1) / ramp2
        inc_dt_ramp2 = int(math.floor(time_ramp2 / dt_nom))
        dt_last_ramp2 = time_ramp2 - inc_dt_ramp2 * dt_nom
        inc_ramp2 = inc_dt_ramp2 if dt_last_ramp2 == 0 else inc_dt_ramp2 + 1
        inc_dt_dwell2 = int(math.floor(time2 / dt_nom))
        dt_last_dwell2 = time2 - inc_dt_dwell2 * dt_nom
        inc_dwell2 = inc_dt_dwell2 if dt_last_dwell2 == 0 else inc_dt_dwell2 + 1
    else:
        dt_last_ramp2 = dt_last_dwell2 = 0
        inc_ramp2 = inc_dwell2 = 0

    if ramp3 > 0:
        time_ramp3 = (temp3 - temp2) / ramp3
        inc_dt_ramp3 = int(math.floor(time_ramp3 / dt_nom))
        dt_last_ramp3 = time_ramp3 - inc_dt_ramp3 * dt_nom
        inc_ramp3 = inc_dt_ramp3 if dt_last_ramp3 == 0 else inc_dt_ramp3 + 1
        inc_dt_dwell3 = int(math.floor(time3 / dt_nom))
        dt_last_dwell3 = time3 - inc_dt_dwell3 * dt_nom
        inc_dwell3 = inc_dt_dwell3 if dt_last_dwell3 == 0 else inc_dt_dwell3 + 1
    else:
        dt_last_ramp3 = dt_last_dwell3 = 0
        inc_ramp3 = inc_dwell3 = 0

    inc = inc_ramp1 + inc_dwell1 + inc_ramp2 + inc_dwell2 + inc_ramp3 + inc_dwell3

    dt = np.zeros(inc + 1)
    dt[1:inc + 1] = dt_nom
    if dt_last_ramp1 > 0:
        dt[inc_ramp1] = dt_last_ramp1
    if dt_last_dwell1 > 0:
        dt[inc_ramp1 + inc_dwell1] = dt_last_dwell1
    if dt_last_ramp2 > 0:
        dt[inc_ramp1 + inc_dwell1 + inc_ramp2] = dt_last_ramp2
    if dt_last_dwell2 > 0:
        dt[inc_ramp1 + inc_dwell1 + inc_ramp2 + inc_dwell2] = dt_last_dwell2
    if dt_last_ramp3 > 0:
        dt[inc_ramp1 + inc_dwell1 + inc_ramp2 + inc_dwell2 + inc_ramp3] = dt_last_ramp3
    if dt_last_dwell3 > 0:
        dt[inc] = dt_last_dwell3

    Time = np.zeros(inc + 1)
    for i in range(1, inc + 1):
        Time[i] = Time[i - 1] + dt[i]

    thermal_profile = np.zeros(inc + 1)
    for i in range(1, inc + 1):
        t = Time[i]
        if ramp2 == 0 and ramp3 == 0:
            thermal_profile[i] = (ramp1 * t + tinit
                                  if t <= (temp1 - tinit) / ramp1 else temp1)
        elif ramp3 == 0:
            if t <= (temp1 - tinit) / ramp1:
                thermal_profile[i] = ramp1 * t + tinit
            elif t <= (temp1 - tinit) / ramp1 + time1:
                thermal_profile[i] = temp1
            elif t <= (temp1 - tinit) / ramp1 + time1 + (temp2 - temp1) / ramp2:
                thermal_profile[i] = ramp2 * (t - ((temp1 - tinit) / ramp1 + time1)) + temp1
            else:
                thermal_profile[i] = temp2
        else:
            if t <= (temp1 - tinit) / ramp1:
                thermal_profile[i] = ramp1 * t + tinit
            elif t <= (temp1 - tinit) / ramp1 + time1:
                thermal_profile[i] = temp1
            elif t <= (temp1 - tinit) / ramp1 + time1 + (temp2 - temp1) / ramp2:
                thermal_profile[i] = ramp2 * (t - ((temp1 - tinit) / ramp1 + time1)) + temp1
            elif t <= (temp1 - tinit) / ramp1 + time1 + (temp2 - temp1) / ramp2 + time2:
                thermal_profile[i] = temp2
            elif t <= ((temp1 - tinit) / ramp1 + time1 + (temp2 - temp1) / ramp2
                       + time2 + (temp3 - temp2) / ramp3):
                thermal_profile[i] = (ramp3 * (t - ((temp1 - tinit) / ramp1 + time1)
                                               - ((temp2 - temp1) / ramp2 + time2)) + temp2)
            else:
                thermal_profile[i] = temp3

    # ---- kinetics parameters ---------------------------------------------
    r = 8.314
    a1, e1, a2, e2 = Cells(2, 11), Cells(3, 11), Cells(4, 11), Cells(5, 11)
    m, n1, n2 = Cells(6, 11), Cells(7, 11), Cells(8, 11)
    # Heat of reaction. Read from Inputs K14, or overridden via the
    # SIMPLECURE_HT environment variable (useful when the sheet's K14 is blank).
    htot = Cells(14, 11)
    _ht_env = os.environ.get("SIMPLECURE_HT")
    if _ht_env not in (None, ""):
        htot = float(_ht_env)
    elif not htot:
        import sys as _sys
        print("WARNING: heat of reaction Cells(14,11)/K14 is empty -> no cure "
              "exotherm (overshoot will be 0). Set SIMPLECURE_HT to override.",
              file=_sys.stderr)
    ad, ed = Cells(9, 11), Cells(10, 11)
    b, w, g = Cells(11, 11), Cells(12, 11), Cells(13, 11)

    a11, e11, a22, e22 = Cells(2, 14), Cells(3, 14), Cells(4, 14), Cells(5, 14)
    m11, n11, m22, n22 = Cells(6, 14), Cells(7, 14), Cells(8, 14), Cells(9, 14)
    d11, aco, act = Cells(10, 14), Cells(11, 14), Cells(12, 14)

    tgo, tgoo, lamda = Cells(18, 11), Cells(19, 11), Cells(20, 11)

    kin_flag = int(Cells(18, 8))
    rom_cp = int(Cells(19, 8))
    rom_k = int(Cells(20, 8))

    def tg_of(alpha):
        return tgo + (tgoo - tgo) * lamda * alpha / (1 - (1 - lamda) * alpha)

    def dadt(tg, temp, conv):
        if kin_flag == 1:
            return dadt1(r, a1, e1, a2, e2, m, n1, n2, ad, ed, b, w, g, tg, temp, conv)
        return dadt2(r, a11, e11, a22, e22, m11, n11, m22, n22, d11, aco, act, temp, conv)

    # ---- constant boundary vectors ---------------------------------------
    # surface convection matrix contribution: h at (nn, nn) only
    hninj_nn = h if convection_top else 0.0

    nj_s3 = np.zeros(nn + 1)
    if convection_top:
        nj_s3[nn] = 1.0

    nj = np.zeros(nn + 1)
    for i in range(1, nn + 1):
        nj[i] = dx / 2 if (i == 1 or i == nn) else dx
    if dx_last > 0:
        nj[nn] = dx_last / 2

    # ---- state arrays -----------------------------------------------------
    cp = np.zeros(ne + 1)
    k_arr = np.zeros(ne + 1)
    t_elm = np.zeros((inc + 1, ne + 1))
    a_elm = np.zeros((inc + 1, ne + 1))
    dadt_elm = np.zeros((inc + 1, ne + 1))
    tg_elm = np.zeros((inc + 1, ne + 1))
    t_node = np.zeros((inc + 1, nn + 1))
    a_node = np.zeros((inc + 1, nn + 1))
    dadt_node = np.zeros((inc + 1, nn + 1))
    tg_node = np.zeros((inc + 1, nn + 1))
    ct_node = np.zeros(nn + 1)

    # ---- time marching ----------------------------------------------------
    for i in range(0, inc + 1):
        if i == 0:
            for j in range(1, nn + 1):
                t_node[i, j] = tinit
                if j <= ne:
                    t_elm[i, j] = tinit
                    a_elm[i, j] = ainit
                    tg_elm[i, j] = tg_of(a_elm[i, j])
                    dadt_elm[i, j] = dadt(tg_elm[i, j], t_elm[i, j], a_elm[i, j])
                ct_node[j] = t_node[i, j]
                if j == 1:
                    dadt_node[i, j] = dadt_elm[i, j]
                if j == nn:
                    dadt_node[i, j] = dadt_elm[i, j - 1]
                if 1 < j < nn:
                    dadt_node[i, j] = (dadt_elm[i, j] + dadt_elm[i, j - 1]) / 2
                a_node[i, j] = ainit
                tg_node[i, j] = tg_of(a_node[i, j])
            continue

        # element cp/k from previous step
        for k1 in range(1, ne + 1):
            cp[k1] = cpc(rom_cp, tg_elm[i - 1, k1], arcp, brcp, drcp, crcp, sigma,
                         afcp, bfcp, wf, t_elm[i - 1, k1])
            k_arr[k1] = kc(rom_k, vf, akr, bkr, ckr, dkr, ekr, fkr, akf, bkf,
                           t_elm[i - 1, k1], a_elm[i - 1, k1])

        # assemble tridiagonal mass (pcpnjni) and stiffness (kdnidnj)
        # NOTE: the VBA 'dxp_last' branches are dead (see module docstring),
        # so no last-element correction is applied here.
        pcpnjni = np.zeros((nn + 1, nn + 1))
        kdnidnj = np.zeros((nn + 1, nn + 1))
        for k1 in range(1, nn + 1):
            for k2 in range(1, nn + 1):
                if abs(k1 - k2) == 1:
                    if k1 > k2:
                        pcpnjni[k1, k2] = p * cp[k1 - 1] * dx / 6
                        kdnidnj[k1, k2] = -k_arr[k1 - 1] / dx
                    else:
                        pcpnjni[k1, k2] = p * cp[k2 - 1] * dx / 6
                        kdnidnj[k1, k2] = -k_arr[k2 - 1] / dx
                elif k1 == k2:
                    if k1 == 1:
                        pcpnjni[k1, k2] = p * cp[k1] * dx / 3
                        kdnidnj[k1, k2] = k_arr[k1] / dx
                    elif k1 == nn:
                        pcpnjni[k1, k2] = p * cp[k1 - 1] * dx / 3
                        kdnidnj[k1, k2] = k_arr[k1 - 1] / dx
                    else:
                        pcpnjni[k1, k2] = p * cp[k1 - 1] * dx / 3 + p * cp[k1] * dx / 3
                        kdnidnj[k1, k2] = k_arr[k1 - 1] / dx + k_arr[k1] / dx

        # surface convection contribution to L
        L = kdnidnj.copy()
        L[nn, nn] += hninj_nn

        mmthdtl = pcpnjni - (1 - theta) * dt[i] * L
        mpthdtl = pcpnjni + theta * dt[i] * L

        # apply Dirichlet rows (identity) to the system matrix
        mpthdtl[1, :] = 0.0
        mpthdtl[1, 1] = 1.0
        if dirichlet_top:
            mpthdtl[nn, :] = 0.0
            mpthdtl[nn, nn] = 1.0

        # load vector cdtf (surface convection + internal cure heat generation)
        cdtf = np.zeros(nn + 1)
        if convection_top:
            sink = too if top_sink_mode == "ambient" else thermal_profile[i]
        else:
            sink = 0.0
        for j in range(1, nn + 1):
            cdtf[j] = (h * sink * dt[i] * nj_s3[j]
                       + dt[i] * nj[j] * htot * pr * (1 - vf) * dadt_node[i - 1, j])

        # RHS: (M-(1-theta)Dt L) t_old + cdtf, with Dirichlet nodes overwritten
        prod = mmthdtl[1:nn + 1, 1:nn + 1] @ ct_node[1:nn + 1]
        cb = np.zeros(nn + 1)
        cb[1:nn + 1] = prod + cdtf[1:nn + 1]
        cb[1] = thermal_profile[i]
        if dirichlet_top:
            cb[nn] = thermal_profile[i]

        ct_node[1:nn + 1] = np.linalg.solve(mpthdtl[1:nn + 1, 1:nn + 1], cb[1:nn + 1])
        t_node[i, 1:nn + 1] = ct_node[1:nn + 1]

        # update elements
        for j in range(1, ne + 1):
            cand = dadt_elm[i - 1, j] * dt[i] + a_elm[i - 1, j]
            a_elm[i, j] = cand if cand <= 1 else 1.0
            t_elm[i, j] = (t_node[i, j] + t_node[i, j + 1]) / 2
            tg_elm[i, j] = tg_of(a_elm[i, j])
            dadt_elm[i, j] = dadt(tg_elm[i, j], t_elm[i, j], a_elm[i, j])

        # update nodal alpha / dadt / Tg
        for j in range(1, nn + 1):
            if j == 1:
                dadt_node[i, j] = dadt_elm[i, j]
                a_node[i, j] = a_elm[i, j]
            if j == nn:
                dadt_node[i, j] = dadt_elm[i, j - 1]
                a_node[i, j] = a_elm[i, j - 1]
            if 1 < j < nn:
                dadt_node[i, j] = (dadt_elm[i, j] + dadt_elm[i, j - 1]) / 2
                a_node[i, j] = (a_elm[i, j] + a_elm[i, j - 1]) / 2
            tg_node[i, j] = tg_of(a_node[i, j])

    # ---- node x-positions (all nn nodes) ---------------------------------
    x_positions = []
    for j in range(1, nn + 1):
        xpos = ((j - 1) * dx) * 1000
        if dx_last > 0 and j == nn:
            xpos = ((j - 2) * dx + dx_last) * 1000
        x_positions.append(xpos)

    time_min = [0.0] + [Time[i] / 60 for i in range(1, inc + 1)]
    temperature = [[t_node[i, j] for j in range(1, nn + 1)] for i in range(0, inc + 1)]
    degree_of_cure = [[a_node[i, j] for j in range(1, nn + 1)] for i in range(0, inc + 1)]
    tg = [[tg_node[i, j] for j in range(1, nn + 1)] for i in range(0, inc + 1)]

    # ---- overshoot / time / location -------------------------------------
    tov2 = 0.0
    timeov2 = 0.0
    inc_over = 0
    over_node = 0
    loc_over = 0.0
    for i in range(1, inc + 1):
        for j in range(1, nn + 1):
            dt_current = t_node[i, j] - thermal_profile[i]
            if dt_current > tov2:
                tov2 = dt_current
                timeov2 = Time[i] / 60
                inc_over = i
                over_node = j
                loc_over = (j - 1) * dx * 1000  # dxp_last branch is dead in VBA

    if inc_over == 0:
        temp_over = 0.0
        time_over = 0.0
    elif inc_over == 1 or inc_over == inc:
        temp_over = tov2
        time_over = timeov2
    else:
        tov1 = t_node[inc_over - 1, over_node] - thermal_profile[inc_over - 1]
        tov3 = t_node[inc_over + 1, over_node] - thermal_profile[inc_over + 1]
        timeov1 = Time[inc_over - 1] / 60
        timeov3 = Time[inc_over + 1] / 60
        alpha = (((tov1 - tov3) * (timeov1 - timeov2)
                  - (tov1 - tov2) * (timeov1 - timeov3))
                 / ((timeov1 ** 2 - timeov3 ** 2) * (timeov1 - timeov2)
                    - (timeov1 - timeov3) * (timeov1 ** 2 - timeov2 ** 2)))
        beta = ((tov1 - tov2) - alpha * (timeov1 ** 2 - timeov2 ** 2)) / (timeov1 - timeov2)
        gamma = tov2 - alpha * timeov2 ** 2 - beta * timeov2
        time_over = -beta / (2 * alpha)
        temp_over = alpha * time_over ** 2 + beta * time_over + gamma

    # ---- cure time --------------------------------------------------------
    threshold = Cells(threshold_row, 2)
    SENTINEL = 8040711.0
    SENTINEL2 = 8045489.0
    ctime2 = 0.0
    ctime_inc = 0
    ctime_node = 0
    for i in range(1, inc + 1):
        amin = SENTINEL
        for j in range(1, nn + 1):
            if tg_node[i, j] < amin:
                amin = tg_node[i, j]
                ctime_node = j
        if amin >= threshold:
            ctime2 = Time[i] / 60
            ctime_inc = i
            break

    if ctime_inc == 0:
        curetime = 0.0
    else:
        ak2 = tg_node[ctime_inc, ctime_node]
        amin2 = SENTINEL2
        for j in range(1, nn + 1):
            if tg_node[ctime_inc - 1, j] < amin2:
                amin2 = tg_node[ctime_inc - 1, j]
        if ctime_inc == 1:
            ak1 = tg_of(ainit)
            ctime1 = 0.0
        else:
            ak1 = amin2
            ctime1 = Time[ctime_inc - 1] / 60
        coeff1 = (ak1 - ak2) / (ctime1 - ctime2)
        coeff2 = ak2 - coeff1 * ctime2
        curetime = (threshold - coeff2) / coeff1

    # ---- maximum temperature ---------------------------------------------
    maxtemp = 0.0
    for i in range(1, inc + 1):
        for j in range(1, nn + 1):
            if t_node[i, j] > maxtemp:
                maxtemp = t_node[i, j]
    if inc_over > 0 and maxtemp < temp_over + thermal_profile[inc_over]:
        maxtemp = temp_over + thermal_profile[inc_over]

    results = {
        "overshoot_C": None if inc_over == 0 else temp_over,
        "overshoot_location_mm": None if inc_over == 0 else loc_over,
        "overshoot_time_min": None if inc_over == 0 else time_over,
        "cure_time_min": None if ctime_inc == 0 else curetime,
        "max_temperature_C": maxtemp,
        "n_increments": inc,
        "n_nodes": nn,
    }

    overshoot_chart = None
    if temp_over > 0:
        overshoot_chart = {
            "time_min": [Time[i] / 60 for i in range(0, inc + 1)],
            "driving_temperature": [tinit] + [thermal_profile[i] for i in range(1, inc + 1)],
            "laminate_temperature": [tinit] + [t_node[i, over_node] for i in range(1, inc + 1)],
            "glass_transition_temperature": [0.0] + [tg_node[i, over_node] for i in range(1, inc + 1)],
            "tg_axis_min": tgo,
            "driving_series_name": series_name_profile,
        }

    return {
        "x_positions": x_positions,
        "time_min": time_min,
        "temperature": temperature,
        "degree_of_cure": degree_of_cure,
        "tg": tg,
        "results": results,
        "overshoot_chart": overshoot_chart,
    }


# ---------------------------------------------------------------------------
# Output helpers (shared)
# ---------------------------------------------------------------------------
def _write_field_csv(path, x_positions, time_min, field):
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["t/x"] + list(x_positions))
        for t, row in zip(time_min, field):
            wr.writerow([t] + list(row))


def write_outputs(out, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    _write_field_csv(os.path.join(output_dir, "Temperature.csv"),
                     out["x_positions"], out["time_min"], out["temperature"])
    _write_field_csv(os.path.join(output_dir, "Degree_of_Cure.csv"),
                     out["x_positions"], out["time_min"], out["degree_of_cure"])
    _write_field_csv(os.path.join(output_dir, "Glass_Transition_Temperature.csv"),
                     out["x_positions"], out["time_min"], out["tg"])

    if out["overshoot_chart"] is not None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            ch = out["overshoot_chart"]
            fig, ax = plt.subplots(figsize=(6, 4.5))
            ax.plot(ch["time_min"], ch["driving_temperature"], label=ch["driving_series_name"])
            ax.plot(ch["time_min"], ch["laminate_temperature"], label="Laminate Temperature")
            ax.plot(ch["time_min"], ch["glass_transition_temperature"],
                    label="Glass Transition Temperature")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel("Temperature [C]")
            ax.set_ylim(bottom=ch["tg_axis_min"])
            ax.legend()
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, "overshoot_chart.png"), dpi=120)
            plt.close(fig)
        except ImportError:
            pass


def print_results(title, results):
    def fmt(v):
        return "N/A" if v is None else f"{v:.6f}"
    print(f"=== {title} -- Main Results ===")
    print(f"  Overshoot [C]           : {fmt(results['overshoot_C'])}")
    print(f"  Overshoot location [mm] : {fmt(results['overshoot_location_mm'])}")
    print(f"  Overshoot time [min]    : {fmt(results['overshoot_time_min'])}")
    print(f"  Cure time [min]         : {fmt(results['cure_time_min'])}")
    print(f"  Maximum Temperature [C] : {fmt(results['max_temperature_C'])}")
    print(f"  (increments={results['n_increments']}, nodes={results['n_nodes']})")
