"""
SimpleCure Autoclave -- Python port of the VBA macro in SimpleCure_Autoclave.xlsm

1-D finite-element (linear elements) transient heat-transfer simulation of the
autoclave cure of a composite laminate on a tool, coupled with resin cure
kinetics and glass-transition-temperature (Tg) evolution.  A theta-method
(Crank-Nicolson-like, theta = 0.999) time integration is used.

This is a faithful translation of:
  * Sheet1.RUN_Click        -> run_simulation()  (main driver)
  * Module1.dadt1           -> dadt1()           (cure kinetics model 1, diffusion)
  * Module7.dadt2           -> dadt2()           (cure kinetics model 2, autocatalytic)
  * Module5.kc              -> kc()              (composite thermal conductivity)
  * Module6.cpc             -> cpc()             (composite specific heat)
  * Module2.inverse / Module3.mmult / Module4.transpose
                            -> replaced by numpy.linalg.solve

Inputs are read from the "Inputs" worksheet of an .xlsm/.xlsx workbook, using the
exact same cell addresses as the VBA (VBA Cells(row, col) with the Inputs sheet
active).  Outputs are written as CSV files and the main results are printed.

Usage:
    python simplecure_autoclave.py [workbook.xlsm] [output_dir]

Defaults: workbook = SimpleCure_Autoclave.xlsm (next to this script),
          output_dir = ./output
"""

import math
import os
import sys
import csv

import numpy as np

try:
    import openpyxl
except ImportError:
    openpyxl = None


# ---------------------------------------------------------------------------
# Module1.dadt1  -- cure kinetics model 1 (with diffusion control)
# ---------------------------------------------------------------------------
def dadt1(r, a1, e1, a2, e2, m, n1, n2, ad, ed, b, w, g, tg, temp, conv):
    """Return da/dt for kinetics model 1 (diffusion-limited)."""
    k1c = a1 * math.exp(-e1 / (r * (temp + 273.15)))
    k2c = a2 * math.exp(-e2 / (r * (temp + 273.15)))

    f = w * (temp - tg) + g

    if f == 0:
        k1 = k1c
        k2 = k2c
        kd = 0.0
    else:
        kd = ad * math.exp(-ed / (r * (temp + 273.15))) * math.exp(-b / f)
        k1 = 0.0 if k1c == 0 else 1.0 / (1.0 / k1c + 1.0 / kd)
        k2 = 0.0 if k2c == 0 else 1.0 / (1.0 / k2c + 1.0 / kd)

    if (1 - conv) > 0.001:
        return k1 * (1 - conv) ** n1 + k2 * conv ** m * (1 - conv) ** n2
    return 0.0


# ---------------------------------------------------------------------------
# Module7.dadt2  -- cure kinetics model 2 (autocatalytic, two-stage)
# ---------------------------------------------------------------------------
def dadt2(r, a11, e11, a22, e22, m11, n11, m22, n22, d11, aco, act, temp, conv):
    """Return da/dt for kinetics model 2."""
    k1 = a11 * math.exp(-e11 / (r * (temp + 273.15)))
    k2 = a22 * math.exp(-e22 / (r * (temp + 273.15)))

    if (1 - conv) > 0.001:
        return (k1 * conv ** m11 * (1 - conv) ** n11
                + (k2 * conv ** m22 * (1 - conv) ** n22)
                / (1 + math.exp(d11 * (conv - (aco + act * (temp + 273.15))))))
    return 0.0


# ---------------------------------------------------------------------------
# Module5.kc  -- composite thermal conductivity (rule of mixtures optional)
# ---------------------------------------------------------------------------
def kc(rom_k, vf, akr, bkr, ckr, dkr, ekr, fkr, akf, bkf, temp, conv):
    kr = (akr * temp * conv + bkr * conv + ckr * temp + dkr
          + ekr * temp * conv ** 2 + fkr * conv ** 2)
    kf = akf * temp + bkf

    if rom_k > 0:
        return (vf * kr * (kf / kr - 1)
                + kr * (0.5 - kf / (2 * kr))
                + kr * (kf / kr - 1)
                * (vf ** 2 - vf + ((kf / kr + 1) ** 2) / ((2 * kf / kr - 2) ** 2)) ** 0.5)
    return kr


# ---------------------------------------------------------------------------
# Module6.cpc  -- composite specific heat (rule of mixtures optional)
# ---------------------------------------------------------------------------
def cpc(rom_cp, tg, arcp, brcp, drcp, crcp, sigma, afcp, bfcp, wf, temp):
    cpr = arcp * temp + brcp + drcp / (1 + math.exp(crcp * (temp - tg - sigma)))
    cpf = afcp * temp + bfcp

    if rom_cp > 0:
        return wf * cpf + (1 - wf) * cpr
    return cpr


# ---------------------------------------------------------------------------
# Input reading -- mimic VBA Cells(row, col) on the "Inputs" sheet
# ---------------------------------------------------------------------------
class Inputs:
    """Grid accessor giving VBA-style 1-based Cells(row, col) into 'Inputs'."""

    def __init__(self, grid):
        # grid[row][col] with 1-based indexing (row/col 0 unused)
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
            raise RuntimeError("openpyxl is required to read the workbook. "
                               "Install with: pip install openpyxl")
        wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
        ws = wb["Inputs"]
        grid = {}
        for r in ws.iter_rows():
            for c in r:
                if c.value is not None:
                    grid.setdefault(c.row, {})[c.column] = c.value
        return cls(grid)


# ---------------------------------------------------------------------------
# Main driver  -- Sheet1.RUN_Click
# ---------------------------------------------------------------------------
def run_simulation(inp):
    """Run the full autoclave cure simulation.

    Parameters
    ----------
    inp : Inputs
        Accessor providing VBA-style Cells(row, col) from the Inputs sheet.

    Returns
    -------
    dict with keys:
        x_positions      : list of laminate depth positions [mm] (part nodes)
        time_min         : list of time stamps [min] (index 0 = initial)
        temperature      : 2D list [time][part-node] temperature [C]
        degree_of_cure   : 2D list [time][part-node] degree of cure
        tg               : 2D list [time][part-node] Tg [C]
        results          : dict of main scalar results
        thermal_profile  : list of autoclave air temperature per increment
    """
    Cells = inp.cell

    # ---- model discretization parameters ----------------------------------
    dx = 0.0005          # nominal element size [m]
    dt_nom = 30          # nominal time increment [s]
    theta = 0.999

    # part elements/nodes
    nep_dx = int(math.floor(Cells(1, 2) / 1000 / dx))
    dxp_last = (Cells(1, 2) / 1000) - (nep_dx * dx)
    nep = nep_dx if dxp_last == 0 else nep_dx + 1
    nnp = nep + 1

    # tool elements/nodes
    net_dx = int(math.floor(Cells(2, 2) / 1000 / dx))
    dxt_last = (Cells(2, 2) / 1000) - (net_dx * dx)
    net = net_dx if dxt_last == 0 else net_dx + 1
    nnt = net + 1

    # totals
    nn = nnp + nnt - 1   # total nodes
    ne = nn - 1          # total elements

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

    vf = Cells(20, 5)               # fibre volume fraction
    p = vf * pf + (1 - vf) * pr     # composite density
    wf = vf * pf / p                # fibre weight fraction

    # ---- tool properties --------------------------------------------------
    ptool = Cells(13, 8)
    cptool = Cells(14, 8)
    ktool = Cells(15, 8)

    # ---- convection -------------------------------------------------------
    htop = Cells(20, 2)
    hbottom = Cells(21, 2)

    # ---- initial conditions -----------------------------------------------
    tinit = Cells(5, 2)
    ainit = Cells(6, 2)

    # ---- thermal profile parameters ---------------------------------------
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

    # ---- increments / time-step calculation -------------------------------
    # first ramp
    time_ramp1 = (temp1 - tinit) / ramp1
    inc_dt_ramp1 = int(math.floor(time_ramp1 / dt_nom))
    dt_last_ramp1 = time_ramp1 - inc_dt_ramp1 * dt_nom
    inc_ramp1 = inc_dt_ramp1 if dt_last_ramp1 == 0 else inc_dt_ramp1 + 1

    # first dwell
    inc_dt_dwell1 = int(math.floor(time1 / dt_nom))
    dt_last_dwell1 = time1 - inc_dt_dwell1 * dt_nom
    inc_dwell1 = inc_dt_dwell1 if dt_last_dwell1 == 0 else inc_dt_dwell1 + 1

    # second ramp + dwell
    if ramp2 > 0:
        time_ramp2 = (temp2 - temp1) / ramp2
        inc_dt_ramp2 = int(math.floor(time_ramp2 / dt_nom))
        dt_last_ramp2 = time_ramp2 - inc_dt_ramp2 * dt_nom
        inc_ramp2 = inc_dt_ramp2 if dt_last_ramp2 == 0 else inc_dt_ramp2 + 1

        inc_dt_dwell2 = int(math.floor(time2 / dt_nom))
        dt_last_dwell2 = time2 - inc_dt_dwell2 * dt_nom
        inc_dwell2 = inc_dt_dwell2 if dt_last_dwell2 == 0 else inc_dt_dwell2 + 1
    else:
        dt_last_ramp2 = 0
        dt_last_dwell2 = 0
        inc_ramp2 = 0
        inc_dwell2 = 0

    # third ramp + dwell
    if ramp3 > 0:
        time_ramp3 = (temp3 - temp2) / ramp3
        inc_dt_ramp3 = int(math.floor(time_ramp3 / dt_nom))
        dt_last_ramp3 = time_ramp3 - inc_dt_ramp3 * dt_nom
        inc_ramp3 = inc_dt_ramp3 if dt_last_ramp3 == 0 else inc_dt_ramp3 + 1

        inc_dt_dwell3 = int(math.floor(time3 / dt_nom))
        dt_last_dwell3 = time3 - inc_dt_dwell3 * dt_nom
        inc_dwell3 = inc_dt_dwell3 if dt_last_dwell3 == 0 else inc_dt_dwell3 + 1
    else:
        dt_last_ramp3 = 0
        dt_last_dwell3 = 0
        inc_ramp3 = 0
        inc_dwell3 = 0

    inc = inc_ramp1 + inc_dwell1 + inc_ramp2 + inc_dwell2 + inc_ramp3 + inc_dwell3

    # time-step array dt(0..inc)
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

    # cumulative time Time(0..inc)
    Time = np.zeros(inc + 1)
    for i in range(1, inc + 1):
        Time[i] = Time[i - 1] + dt[i]

    # thermal profile (autoclave air temperature) per increment 1..inc
    thermal_profile = np.zeros(inc + 1)
    for i in range(1, inc + 1):
        t = Time[i]
        if ramp2 == 0 and ramp3 == 0:            # 1 dwell cycle
            if t <= (temp1 - tinit) / ramp1:
                thermal_profile[i] = ramp1 * t + tinit
            else:
                thermal_profile[i] = temp1
        elif ramp3 == 0:                          # 2 dwell cycle
            if t <= (temp1 - tinit) / ramp1:
                thermal_profile[i] = ramp1 * t + tinit
            elif t <= (temp1 - tinit) / ramp1 + time1:
                thermal_profile[i] = temp1
            elif t <= (temp1 - tinit) / ramp1 + time1 + (temp2 - temp1) / ramp2:
                thermal_profile[i] = ramp2 * (t - ((temp1 - tinit) / ramp1 + time1)) + temp1
            else:
                thermal_profile[i] = temp2
        else:                                     # 3 dwell cycle
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

    # ---- kinetic parameters (model 1) -------------------------------------
    r = 8.314
    a1, e1 = Cells(2, 11), Cells(3, 11)
    a2, e2 = Cells(4, 11), Cells(5, 11)
    m = Cells(6, 11)
    n1, n2 = Cells(7, 11), Cells(8, 11)
    htot = Cells(14, 11)
    ad, ed = Cells(9, 11), Cells(10, 11)
    b, w, g = Cells(11, 11), Cells(12, 11), Cells(13, 11)

    # ---- kinetic parameters (model 2) -------------------------------------
    a11, e11 = Cells(2, 14), Cells(3, 14)
    a22, e22 = Cells(4, 14), Cells(5, 14)
    m11, n11 = Cells(6, 14), Cells(7, 14)
    m22, n22 = Cells(8, 14), Cells(9, 14)
    d11 = Cells(10, 14)
    aco, act = Cells(11, 14), Cells(12, 14)

    # ---- glass transition parameters --------------------------------------
    tgo = Cells(18, 11)
    tgoo = Cells(19, 11)
    lamda = Cells(20, 11)

    # ---- flags ------------------------------------------------------------
    kin_flag = int(Cells(18, 8))
    rom_cp = int(Cells(19, 8))
    rom_k = int(Cells(20, 8))

    def tg_of(alpha):
        """DiBenedetto Tg as a function of degree of cure alpha."""
        return tgo + (tgoo - tgo) * lamda * alpha / (1 - (1 - lamda) * alpha)

    def dadt(tg, temp, conv):
        if kin_flag == 1:
            return dadt1(r, a1, e1, a2, e2, m, n1, n2, ad, ed, b, w, g, tg, temp, conv)
        return dadt2(r, a11, e11, a22, e22, m11, n11, m22, n22, d11, aco, act, temp, conv)

    # ---- constant boundary matrices (1-based, index 1..nn) ----------------
    hninj_s3 = np.zeros((nn + 1, nn + 1))
    hninj_s3[nn, nn] = htop        # convection on top (part surface)
    hninj_s3[1, 1] = hbottom       # convection on bottom (tool surface)

    nj_s3 = np.zeros(nn + 1)
    nj_s3[1] = 1.0
    nj_s3[nn] = 1.0

    nj = np.zeros(nn + 1)
    for i in range(1, nn + 1):
        # part
        if i == nn or i == nnt:
            nj[i] = dx / 2
        else:
            nj[i] = dx
        # tool nodes zeroed (heat generation only in part)
        if i < nnt:
            nj[i] = 0.0
    if dxp_last > 0:
        nj[nn] = dxp_last / 2

    htc = np.zeros(nn + 1)
    htc[1] = hbottom
    htc[nn] = htop

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

    # ct_node holds the "current" nodal temperature vector carried between steps
    ct_node = np.zeros(nn + 1)

    # ---- time-space marching ---------------------------------------------
    for i in range(0, inc + 1):
        if i == 0:
            # initial conditions
            for j in range(1, nn + 1):
                t_node[i, j] = tinit
                if j <= ne:
                    t_elm[i, j] = tinit
                    if j <= nnt - 1:
                        a_elm[i, j] = 0.0      # tool elements
                        dadt_elm[i, j] = 0.0
                        tg_elm[i, j] = 0.0
                    else:
                        a_elm[i, j] = ainit    # part elements
                        tg_elm[i, j] = tg_of(a_elm[i, j])
                        dadt_elm[i, j] = dadt(tg_elm[i, j], t_elm[i, j], a_elm[i, j])

                ct_node[j] = t_node[i, j]

                if j == nnt:
                    dadt_node[i, j] = dadt_elm[i, j]
                if j == nn:
                    dadt_node[i, j] = dadt_elm[i, j - 1]
                if nnt < j < nn:
                    dadt_node[i, j] = (dadt_elm[i, j] + dadt_elm[i, j - 1]) / 2

                if j < nnt:
                    a_node[i, j] = 0.0
                    dadt_node[i, j] = 0.0
                    tg_node[i, j] = 0.0
                else:
                    a_node[i, j] = ainit
                    tg_node[i, j] = tg_of(a_node[i, j])
            continue

        # ----- i > 0 : build element cp/k from previous step -----
        for k1 in range(1, ne + 1):
            cp[k1] = cpc(rom_cp, tg_elm[i - 1, k1], arcp, brcp, drcp, crcp, sigma,
                         afcp, bfcp, wf, t_elm[i - 1, k1])
            k_arr[k1] = kc(rom_k, vf, akr, bkr, ckr, dkr, ekr, fkr, akf, bkf,
                           t_elm[i - 1, k1], a_elm[i - 1, k1])

        # ----- assemble global mass (pcpnjni) and stiffness (kdnidnj) -----
        pcpnjni = np.zeros((nn + 1, nn + 1))
        kdnidnj = np.zeros((nn + 1, nn + 1))

        for k1 in range(1, nn + 1):
            for k2 in range(1, nn + 1):
                # ---- tool block ----
                if k1 <= nnt and k2 <= nnt:
                    if abs(k1 - k2) == 1:
                        pcpnjni[k1, k2] = ptool * cptool * dx / 6
                        kdnidnj[k1, k2] = -ktool / dx
                        if dxt_last > 0:
                            if k1 == nnt and k2 == nnt - 1:
                                pcpnjni[k1, k2] = ptool * cptool * dxt_last / 6
                                kdnidnj[k1, k2] = -ktool / dxt_last
                            if k2 == nnt and k1 == nnt - 1:
                                pcpnjni[k1, k2] = ptool * cptool * dxt_last / 6
                                kdnidnj[k1, k2] = -ktool / dxt_last
                    elif k1 == k2:
                        if k1 == 1 or k2 == nn:
                            pcpnjni[k1, k2] = ptool * cptool * dx / 3
                            kdnidnj[k1, k2] = ktool / dx
                        else:
                            pcpnjni[k1, k2] = ptool * cptool * dx / 1.5
                            kdnidnj[k1, k2] = 2 * ktool / dx
                        if dxt_last > 0:
                            if k1 == nnt - 1:
                                pcpnjni[k1, k2] = (ptool * cptool * dx / 3
                                                   + ptool * cptool * dxt_last / 3)
                                kdnidnj[k1, k2] = ktool / dx + ktool / dxt_last
                    else:
                        pcpnjni[k1, k2] = 0.0
                        kdnidnj[k1, k2] = 0.0

                # ---- part block ----
                if k1 >= nnt and k2 >= nnt:
                    if abs(k1 - k2) == 1:
                        if k1 > k2:
                            pcpnjni[k1, k2] = p * cp[k1 - 1] * dx / 6
                            kdnidnj[k1, k2] = -k_arr[k1 - 1] / dx
                        else:
                            pcpnjni[k1, k2] = p * cp[k2 - 1] * dx / 6
                            kdnidnj[k1, k2] = -k_arr[k2 - 1] / dx
                        if dxp_last > 0:
                            if k1 == nn and k2 == nn - 1:
                                pcpnjni[k1, k2] = p * cp[k1 - 1] * dxp_last / 6
                                kdnidnj[k1, k2] = -k_arr[k1 - 1] / dxp_last
                            if k2 == nn and k1 == nn - 1:
                                pcpnjni[k1, k2] = p * cp[k2 - 1] * dxp_last / 6
                                kdnidnj[k1, k2] = -k_arr[k2 - 1] / dxp_last
                    elif k1 == k2:
                        if k1 == 1:
                            pcpnjni[k1, k2] = p * cp[k1] * dx / 3
                            kdnidnj[k1, k2] = k_arr[k1] / dx
                        elif k1 == nn:
                            if dxp_last > 0:
                                pcpnjni[k1, k2] = p * cp[k1 - 1] * dxp_last / 3
                                kdnidnj[k1, k2] = k_arr[k1 - 1] / dxp_last
                            else:
                                pcpnjni[k1, k2] = p * cp[k1 - 1] * dx / 3
                                kdnidnj[k1, k2] = k_arr[k1 - 1] / dx
                        else:
                            pcpnjni[k1, k2] = (p * cp[k1 - 1] * dx / 3
                                               + p * cp[k1] * dx / 3)
                            kdnidnj[k1, k2] = k_arr[k1 - 1] / dx + k_arr[k1] / dx
                            if dxp_last > 0:
                                if k1 == nn - 1:
                                    pcpnjni[k1, k2] = (p * cp[k1 - 1] * dx / 3
                                                       + p * cp[k1] * dxp_last / 3)
                                    kdnidnj[k1, k2] = (k_arr[k1 - 1] / dx
                                                       + k_arr[k1] / dxp_last)
                    else:
                        pcpnjni[k1, k2] = 0.0
                        kdnidnj[k1, k2] = 0.0

                # ---- shared tool/part node ----
                if k1 == nnt and k2 == nnt:
                    if dxt_last > 0:
                        pcpnjni[k1, k2] = (p * cp[k1] * dx / 3
                                           + ptool * cptool * dxt_last / 3)
                        kdnidnj[k1, k2] = k_arr[k1] / dx + ktool / dxt_last
                    else:
                        pcpnjni[k1, k2] = (p * cp[k1] * dx / 3
                                           + ptool * cptool * dx / 3)
                        kdnidnj[k1, k2] = k_arr[k1] / dx + ktool / dx

        # M + theta*Dt*L   and   M - (1-theta)*Dt*L
        L = kdnidnj + hninj_s3
        mpthdtl = pcpnjni + theta * dt[i] * L
        mmthdtl = pcpnjni - (1 - theta) * dt[i] * L

        # load / source vector
        cdtf = np.zeros(nn + 1)
        for j in range(1, nn + 1):
            cdtf[j] = (htc[j] * thermal_profile[i] * dt[i] * nj_s3[j]
                       + dt[i] * nj[j] * htot * pr * (1 - vf) * dadt_node[i - 1, j])

        # cb = (M-(1-theta)Dt L) * t_old + cdtf
        ccb = mmthdtl[1:nn + 1, 1:nn + 1] @ ct_node[1:nn + 1]
        cb = np.zeros(nn + 1)
        cb[1:nn + 1] = ccb + cdtf[1:nn + 1]

        # solve (M+theta Dt L) t_new = cb   (equivalent to inverse * cb)
        ct_node[1:nn + 1] = np.linalg.solve(mpthdtl[1:nn + 1, 1:nn + 1], cb[1:nn + 1])

        t_node[i, 1:nn + 1] = ct_node[1:nn + 1]

        # ----- update element quantities -----
        for j in range(1, ne + 1):
            if j <= nnt - 1:
                a_elm[i, j] = 0.0
                tg_elm[i, j] = 0.0
            else:
                cand = dadt_elm[i - 1, j] * dt[i] + a_elm[i - 1, j]
                a_elm[i, j] = cand if cand <= 1 else 1.0
                tg_elm[i, j] = tg_of(a_elm[i, j])
            t_elm[i, j] = (t_node[i, j] + t_node[i, j + 1]) / 2
            if j <= nnt - 1:
                dadt_elm[i, j] = 0.0
            else:
                dadt_elm[i, j] = dadt(tg_elm[i, j], t_elm[i, j], a_elm[i, j])

        # ----- update nodal alpha / dadt / Tg -----
        for j in range(1, nn + 1):
            if j == nnt:
                dadt_node[i, j] = dadt_elm[i, j]
                a_node[i, j] = a_elm[i, j]
            if j == nn:
                dadt_node[i, j] = dadt_elm[i, j - 1]
                a_node[i, j] = a_elm[i, j - 1]
            if nnt < j < nn:
                dadt_node[i, j] = (dadt_elm[i, j] + dadt_elm[i, j - 1]) / 2
                a_node[i, j] = (a_elm[i, j] + a_elm[i, j - 1]) / 2

            if j < nnt:
                dadt_node[i, j] = 0.0
                a_node[i, j] = 0.0
                tg_node[i, j] = 0.0
            else:
                tg_node[i, j] = tg_of(a_node[i, j])

    # -----------------------------------------------------------------------
    # Post-processing:  laminate depth positions (part nodes nnt..nn)
    # -----------------------------------------------------------------------
    x_positions = []
    for j in range(1, nnp + 1):
        xpos = ((j - 1) * dx) * 1000
        if dxp_last > 0 and j == nnp:
            xpos = ((j - 2) * dx + dxp_last) * 1000
        x_positions.append(xpos)

    # time stamps (index 0 = initial state, then increments 1..inc)
    time_min = [0.0] + [Time[i] / 60 for i in range(1, inc + 1)]

    temperature = []
    degree_of_cure = []
    tg = []
    for i in range(0, inc + 1):
        temperature.append([t_node[i, j] for j in range(nnt, nn + 1)])
        degree_of_cure.append([a_node[i, j] for j in range(nnt, nn + 1)])
        tg.append([tg_node[i, j] for j in range(nnt, nn + 1)])

    # -----------------------------------------------------------------------
    # temperature overshoot / time / location  (quadratic peak interpolation)
    # -----------------------------------------------------------------------
    tov2 = 0.0
    timeov2 = 0.0
    inc_over = 0
    over_node = 0
    loc_over = 0.0
    for i in range(1, inc + 1):
        for j in range(nnt, nn + 1):
            dt_current = t_node[i, j] - thermal_profile[i]
            if dt_current > tov2:
                tov2 = dt_current
                timeov2 = Time[i] / 60
                inc_over = i
                over_node = j
                if dxp_last > 0 and j == nn:
                    loc_over = (j - 2) * dx * 1000 + dxp_last - Cells(2, 2)
                else:
                    loc_over = (j - 1) * dx * 1000 - Cells(2, 2)

    if inc_over == 0:
        temp_over = 0.0
        time_over = 0.0
    else:
        if inc_over == 1 or inc_over == inc:
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

    # -----------------------------------------------------------------------
    # cure time  (time when minimum part Tg reaches the threshold)
    # -----------------------------------------------------------------------
    threshold = Cells(23, 2)
    SENTINEL = 8040711.0
    SENTINEL2 = 8045489.0

    ctime2 = 0.0
    ctime_inc = 0
    ctime_node = 0
    for i in range(1, inc + 1):
        amin = SENTINEL
        for j in range(nnt, nn + 1):
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
        for j in range(nnt, nn + 1):
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

    # -----------------------------------------------------------------------
    # maximum temperature
    # -----------------------------------------------------------------------
    maxtemp = 0.0
    for i in range(1, inc + 1):
        for j in range(nnt, nn + 1):
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
        "n_part_nodes": nnp,
        "n_tool_nodes": nnt,
    }

    # overshoot chart series (matches the VBA chart at the overshoot node)
    overshoot_chart = None
    if temp_over > 0:
        time_list = [Time[i] / 60 for i in range(0, inc + 1)]
        laminate = [tinit] + [t_node[i, over_node] for i in range(1, inc + 1)]
        air = [tinit] + [thermal_profile[i] for i in range(1, inc + 1)]
        tg_series = [0.0] + [tg_node[i, over_node] for i in range(1, inc + 1)]
        overshoot_chart = {
            "time_min": time_list,
            "autoclave_temperature": air,
            "laminate_temperature": laminate,
            "glass_transition_temperature": tg_series,
            "tg_axis_min": tgo,
        }

    return {
        "x_positions": x_positions,
        "time_min": time_min,
        "temperature": temperature,
        "degree_of_cure": degree_of_cure,
        "tg": tg,
        "results": results,
        "thermal_profile": [0.0] + list(thermal_profile[1:inc + 1]),
        "overshoot_chart": overshoot_chart,
    }


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
def _write_field_csv(path, x_positions, time_min, field):
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["t/x"] + [x for x in x_positions])
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

    # optional overshoot chart PNG
    if out["overshoot_chart"] is not None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            ch = out["overshoot_chart"]
            fig, ax = plt.subplots(figsize=(6, 4.5))
            ax.plot(ch["time_min"], ch["autoclave_temperature"], label="Autoclave Temperature")
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
            pass  # matplotlib optional


def print_results(results):
    def fmt(v):
        return "N/A" if v is None else f"{v:.6f}"
    print("=== Main Results ===")
    print(f"  Overshoot [C]           : {fmt(results['overshoot_C'])}")
    print(f"  Overshoot location [mm] : {fmt(results['overshoot_location_mm'])}")
    print(f"  Overshoot time [min]    : {fmt(results['overshoot_time_min'])}")
    print(f"  Cure time [min]         : {fmt(results['cure_time_min'])}")
    print(f"  Maximum Temperature [C] : {fmt(results['max_temperature_C'])}")
    print(f"  (increments={results['n_increments']}, nodes={results['n_nodes']}, "
          f"part nodes={results['n_part_nodes']}, tool nodes={results['n_tool_nodes']})")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    workbook = argv[1] if len(argv) > 1 else os.path.join(here, "SimpleCure_Autoclave.xlsm")
    output_dir = argv[2] if len(argv) > 2 else os.path.join(here, "output")

    inp = Inputs.from_workbook(workbook)
    out = run_simulation(inp)
    print_results(out["results"])
    write_outputs(out, output_dir)
    print(f"\nField CSVs and chart written to: {output_dir}")


if __name__ == "__main__":
    main(sys.argv)
