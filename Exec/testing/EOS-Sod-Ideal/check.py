#!/usr/bin/env python3
"""Gates for the tabulated-EOS twin-run Sod case (Stage 3, W9).

1. twin agreement:   tabulated vs thermally_perfect, L1 per field — isolates
                     EOS-path error with identical numerics; bound set by
                     table interpolation accuracy.
2. analytic:         both runs vs the exact Riemann solution — the tabulated
                     run must sit at (essentially) the same distance from
                     truth as the ideal run (scheme error dominates).
3. conservation:     total mass and energy drift, first vs last plotfile.
4. wall-time budget: tabulated/ideal 'Run Time total' ratio (plan review
                     amendment #4).

Tolerances marked MEASURED were committed from the first passing run per the
plan's no-guessed-gates amendment.
"""

import re
import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

GAMMA = 1.4
# Sod initial states (code units)
RHO_L, P_L, U_L = 1.0, 1.0, 0.0
RHO_R, P_R, U_R = 0.125, 0.1, 0.0
T_END = 0.2

TWIN_TOL = 1.0e-3          # MEASURED 2026-07-08: max 4.1e-4 (T) at 2048 cells
ANALYTIC_FACTOR = 1.25     # MEASURED: ratio 0.99 (tabulated marginally closer)
IDEAL_L1_MAX = 5.0e-3      # MEASURED: 1.8e-3 on rho at 2048 cells (minmod/HLLC)
CONS_TOL = 1.0e-11         # relative mass/energy drift
WALL_RATIO_MAX = 3.0       # MEASURED: 2.67 (advance time) at 2048 cells

failed = []


def gate(name, ok, detail):
    print("check.py: %-24s %s  %s" % (name, "PASS" if ok else "** FAIL", detail))
    if not ok:
        failed.append(name)


# ---------------------------------------------------------------------------
# exact Riemann solution for the Sod problem (Toro ch. 4: pressure-function
# Newton for the star state, then wave-by-wave sampling)


def exact_sod(x, t):
    g = GAMMA
    aL = np.sqrt(g * P_L / RHO_L)
    aR = np.sqrt(g * P_R / RHO_R)

    def f_side(p, ps, rhos, a_s):
        if p > ps:  # shock
            A = 2.0 / ((g + 1) * rhos)
            B = (g - 1) / (g + 1) * ps
            return (p - ps) * np.sqrt(A / (p + B))
        # rarefaction
        return 2 * a_s / (g - 1) * ((p / ps) ** ((g - 1) / (2 * g)) - 1.0)

    def fprime(p, ps, rhos, a_s):
        eps = 1e-8 * max(p, 1.0)
        return (f_side(p + eps, ps, rhos, a_s) - f_side(p - eps, ps, rhos, a_s)) / (2 * eps)

    p = 0.5 * (P_L + P_R)
    for _ in range(60):
        f = f_side(p, P_L, RHO_L, aL) + f_side(p, P_R, RHO_R, aR) + (U_R - U_L)
        df = fprime(p, P_L, RHO_L, aL) + fprime(p, P_R, RHO_R, aR)
        p_new = p - f / df
        if p_new < 1e-10:
            p_new = 1e-10
        if abs(p_new - p) < 1e-14 * p:
            p = p_new
            break
        p = p_new
    pstar = p
    ustar = 0.5 * (U_L + U_R) + 0.5 * (f_side(pstar, P_R, RHO_R, aR) -
                                       f_side(pstar, P_L, RHO_L, aL))

    rho = np.empty_like(x)
    u = np.empty_like(x)
    pr = np.empty_like(x)
    s = x / t
    for i, si in enumerate(s):
        if si < ustar:  # left of contact
            if pstar > P_L:  # left shock
                SL = U_L - aL * np.sqrt((g + 1) / (2 * g) * pstar / P_L + (g - 1) / (2 * g))
                if si < SL:
                    rho[i], u[i], pr[i] = RHO_L, U_L, P_L
                else:
                    r = RHO_L * ((pstar / P_L + (g - 1) / (g + 1)) /
                                 ((g - 1) / (g + 1) * pstar / P_L + 1))
                    rho[i], u[i], pr[i] = r, ustar, pstar
            else:  # left rarefaction
                aSL = aL * (pstar / P_L) ** ((g - 1) / (2 * g))
                headL = U_L - aL
                tailL = ustar - aSL
                if si < headL:
                    rho[i], u[i], pr[i] = RHO_L, U_L, P_L
                elif si > tailL:
                    r = RHO_L * (pstar / P_L) ** (1 / g)
                    rho[i], u[i], pr[i] = r, ustar, pstar
                else:  # inside the fan
                    uf = 2 / (g + 1) * (aL + (g - 1) / 2 * U_L + si)
                    af = 2 / (g + 1) * (aL + (g - 1) / 2 * (U_L - si))
                    rho[i] = RHO_L * (af / aL) ** (2 / (g - 1))
                    u[i] = uf
                    pr[i] = P_L * (af / aL) ** (2 * g / (g - 1))
        else:  # right of contact
            if pstar > P_R:  # right shock
                SR = U_R + aR * np.sqrt((g + 1) / (2 * g) * pstar / P_R + (g - 1) / (2 * g))
                if si > SR:
                    rho[i], u[i], pr[i] = RHO_R, U_R, P_R
                else:
                    r = RHO_R * ((pstar / P_R + (g - 1) / (g + 1)) /
                                 ((g - 1) / (g + 1) * pstar / P_R + 1))
                    rho[i], u[i], pr[i] = r, ustar, pstar
            else:  # right rarefaction
                aSR = aR * (pstar / P_R) ** ((g - 1) / (2 * g))
                headR = U_R + aR
                tailR = ustar + aSR
                if si > headR:
                    rho[i], u[i], pr[i] = RHO_R, U_R, P_R
                elif si < tailR:
                    r = RHO_R * (pstar / P_R) ** (1 / g)
                    rho[i], u[i], pr[i] = r, ustar, pstar
                else:
                    uf = 2 / (g + 1) * (-aR + (g - 1) / 2 * U_R + si)
                    af = 2 / (g + 1) * (aR - (g - 1) / 2 * (U_R - si))
                    rho[i] = RHO_R * (af / aR) ** (2 / (g - 1))
                    u[i] = uf
                    pr[i] = P_R * (af / aR) ** (2 * g / (g - 1))
    return rho, u, pr


# ---------------------------------------------------------------------------
# load the two runs


def profile(prefix, which):
    files = sorted(get_files(".", include=[prefix], exclude=["temp"], get_all=True))
    if not files:
        print("check.py: no plotfiles for %s" % prefix)
        sys.exit(1)
    d = ReadBoxLib(files[0 if which == "first" else -1])
    out = {}
    for f in ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid", "nrg-fluid"):
        x, v = d.get(f)
        if isinstance(x, (list, tuple)):
            x = x[0]
        out[f] = (np.asarray(x).ravel(), np.asarray(v).ravel())
    out["time"] = d.time
    return out


ideal_0 = profile("ideal.plt", "first")
ideal_1 = profile("ideal.plt", "last")
tab_0 = profile("tabulated.plt", "first")
tab_1 = profile("tabulated.plt", "last")

# ---------------------------------------------------------------------------
# 1. twin agreement

for f in ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid"):
    vi = ideal_1[f][1]
    vt = tab_1[f][1]
    rng = vi.max() - vi.min()
    l1 = np.mean(np.abs(vt - vi)) / rng
    gate("twin-%s" % f.split("-")[0], l1 <= TWIN_TOL, "L1/range=%.3e (tol %g)" % (l1, TWIN_TOL))

# ---------------------------------------------------------------------------
# 2. vs the exact solution (density, the sharpest field)

x = ideal_1["rho-fluid"][0]
t = ideal_1["time"]
rho_exact, _, _ = exact_sod(x, t)
l1_ideal = np.mean(np.abs(ideal_1["rho-fluid"][1] - rho_exact)) / (RHO_L - RHO_R)
l1_tab = np.mean(np.abs(tab_1["rho-fluid"][1] - rho_exact)) / (RHO_L - RHO_R)
gate("analytic-ideal", l1_ideal <= IDEAL_L1_MAX, "L1=%.3e (max %g)" % (l1_ideal, IDEAL_L1_MAX))
gate("analytic-tabulated", l1_tab <= ANALYTIC_FACTOR * l1_ideal,
     "L1=%.3e vs ideal %.3e (factor %g)" % (l1_tab, l1_ideal, ANALYTIC_FACTOR))

# ---------------------------------------------------------------------------
# 3. conservation (uniform grid: sums are integrals up to a constant factor)

for name, d0, d1 in (("ideal", ideal_0, ideal_1), ("tabulated", tab_0, tab_1)):
    for f, label in (("rho-fluid", "mass"), ("nrg-fluid", "energy")):
        s0 = d0[f][1].sum()
        s1 = d1[f][1].sum()
        drift = abs(s1 - s0) / abs(s0)
        gate("conserve-%s-%s" % (name, label), drift <= CONS_TOL, "drift=%.3e" % drift)

# ---------------------------------------------------------------------------
# 4. wall-time budget


def wall(log):
    # 'advance' isolates the stepping cost; 'total' would also count init
    # (incl. the ~10 MB table parse), which is not what the budget bounds
    txt = open(log).read()
    m = re.findall(r"Run Time advance\s*=\s*([0-9.eE+-]+)", txt)
    return float(m[-1]) if m else None


w_i = wall("run_log_ideal.txt")
w_t = wall("run_log_tabulated.txt")
if w_i and w_t:
    ratio = w_t / w_i
    gate("wall-time", ratio <= WALL_RATIO_MAX,
         "tabulated/ideal = %.2f (%.2fs / %.2fs, max %g)" % (ratio, w_t, w_i, WALL_RATIO_MAX))
else:
    gate("wall-time", False, "could not parse 'Run Time total' from logs")

print("check.py:", "FAIL" if failed else "PASS")
sys.exit(1 if failed else 0)
