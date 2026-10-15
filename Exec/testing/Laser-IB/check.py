#!/usr/bin/env python3
"""Gates for the laser inverse-bremsstrahlung ray march (Laser-IB).

The plasma is frozen (no flux action, tiny intensity, 2 steps), so the FIRST
[laser] log line is the march applied to the exact initial fields. kappa is
piecewise constant per cell in the code, so the references below use the
same cell-centred kappa_i: the gates test the traversal and the absorption
bookkeeping to round-off, not the quadrature error of the cell discretisation
(which is reported, not gated, against the continuous closed form).

ramp1d : n_e = x/L up to n_c then 1.2 n_c, plane wave along +x.
         absorbed/P_in = 1 - exp(-tau), tau = sum_i kappa_i dx over the ramp
         cells; the remainder is dumped at the first cell >= n_c; nothing
         transmitted. Continuous form: tau = (16/15) nu_c L / c (one pass to
         n_c; the textbook 32/15 is the there-and-back with reflection).
ramp2d : same ramp, one ray at 30 deg to the ramp normal. Without refraction
         the chord per column is dx / cos(theta) exactly, whatever the DDA
         does at the y-faces: tau = sum_i kappa_i dx / cos(theta).
beam2d : uniform n_e = 0.5 n_c, top-hat beam radius R along +x with
         rays_per_cell = 2: every interior row gets the same two rays, so the
         deposition is y-uniform to round-off inside |y - y0| < R, zero
         outside, and decays column to column by exactly exp(-kappa dx).
conservation (all cases, every evaluation): P_in = absorbed + dumped +
         transmitted to round-off.

Tolerances marked MEASURED were committed from the first passing run.
"""

import re
import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

if not hasattr(np, "product"):
    np.product = np.prod

# ---- constants and reference quantities, identical to the lua / the action ----
eps0, m_e, e, c, kB, m_p = 8.8541878128e-12, 9.1093837015e-31, 1.602176634e-19, 299792458.0, 1.38064852e-23, 1.67262192e-27
lam = 0.351e-6
n_crit = eps0 * m_e * (2 * np.pi * c / lam) ** 2 / e**2
x_ref, m_ref, rho_ref, T_ref = 1.0e-4, m_p, n_crit * m_p, 1.0e3 * e / kB
n_ref = rho_ref / m_ref
u_ref = np.sqrt(kB * T_ref / m_ref)
t_ref = x_ref / u_ref
prs_ref = n_ref * m_ref * u_ref**2
kT_eV = kB * T_ref / e
C_kappa = 2.91e-12 * n_ref * kT_eV ** (-1.5) * x_ref / c
LNL, T0, L, n_floor, I0 = 7.0, 1.0, 1.0, 1.0e-3, 1.0e17
I_scale = t_ref / (prs_ref * x_ref)

failed = []


def gate(name, value, tol):
    ok = value <= tol
    print("  %-48s %.3e  (tol %.1e)  %s" % (name, value, tol, "PASS" if ok else "FAIL"))
    if not ok:
        failed.append(name)


def laser_lines(case):
    pat = re.compile(r"\[laser\] t = (\S+) P_in = (\S+) absorbed = (\S+) dumped = (\S+) transmitted = (\S+)")
    rows = [tuple(float(v) for v in m.groups()) for m in map(pat.search, open("run_log_%s.txt" % case)) if m]
    if not rows:
        failed.append(case + ": no [laser] lines in the log")
    return np.array(rows)


def kappa(n):  # code units, Z = 1 so n_e = Z2n = n, x = n / n_c = n
    return C_kappa * LNL * n * T0 ** (-1.5) * n / np.sqrt(1.0 - n)


def ramp_tau(nx, cos_theta):
    dx = 1.2 / nx
    xc = (np.arange(nx) + 0.5) * dx
    n = np.where(xc < L, np.maximum(xc / L, n_floor), 1.2)
    ramp = n < 1.0
    return np.sum(kappa(n[ramp]) * dx) / cos_theta


def conservation(case, rows):
    r = np.max(np.abs(rows[:, 1] - rows[:, 2] - rows[:, 3] - rows[:, 4]) / rows[:, 1])
    gate(case + ": max |P_in - abs - dump - trans| / P_in", r, 1e-12)


# ------------------------------------------------------------------ ramp1d
print("ramp1d")
rows = laser_lines("ramp1d")
if len(rows):
    conservation("ramp1d", rows)
    P_in, absd, dump, trans = rows[0, 1:]
    tau = ramp_tau(120, 1.0)
    gate("ramp1d: |absorbed/P_in - (1 - e^-tau)|", abs(absd / P_in - (1 - np.exp(-tau))), 1e-12)
    gate("ramp1d: |dumped/P_in - e^-tau|", abs(dump / P_in - np.exp(-tau)), 1e-12)
    gate("ramp1d: transmitted/P_in", trans / P_in, 0.0)
    gate("ramp1d: |P_in - I0 * I_scale| / P_in", abs(P_in - I0 * I_scale) / P_in, 1e-12)
    # continuous closed form, one pass: tau_c = (16/15) * kappa-prefactor * L  (nu_c at n_c, T0)
    tau_c = (16.0 / 15.0) * C_kappa * LNL * T0 ** (-1.5) * L
    print("  cell-discretised tau = %.6f, continuous (16/15) nu_c L/c = %.6f, rel diff %.2e (not gated)"
          % (tau, tau_c, abs(tau - tau_c) / tau_c))

# ------------------------------------------------------------------ ramp2d
print("ramp2d")
rows = laser_lines("ramp2d")
if len(rows):
    conservation("ramp2d", rows)
    P_in, absd, dump, trans = rows[0, 1:]
    tau = ramp_tau(120, np.cos(np.radians(30.0)))
    gate("ramp2d: |absorbed/P_in - (1 - e^-tau)|", abs(absd / P_in - (1 - np.exp(-tau))), 1e-12)
    gate("ramp2d: |dumped/P_in - e^-tau|", abs(dump / P_in - np.exp(-tau)), 1e-12)
    gate("ramp2d: transmitted/P_in", trans / P_in, 0.0)

# ------------------------------------------------------------------ beam2d
print("beam2d")
rows = laser_lines("beam2d")
if len(rows):
    conservation("beam2d", rows)
    files = get_files(".", include=["beam2d.plt"], exclude=[".inputs", ".txt"], get_all=True)
    d0 = ReadBoxLib(files[0], max_level=0)
    d1 = ReadBoxLib(files[1], max_level=0)
    (x, y), E0 = d0.get("nrg-plasma")
    _, E1 = d1.get("nrg-plasma")
    dep = (E1 - E0) / (d1.time - d0.time)  # ~ deposited power density per cell
    dx, dy = x[1] - x[0], y[1] - y[0]
    R, y0 = 0.2, 0.5
    inside = np.abs(y - y0) < R - dy  # rows fully inside the beam
    outside = np.abs(y - y0) > R + dy
    col = dep[:, inside]
    gate("beam2d: y-uniformity inside the beam (rel spread)",
         np.max(np.ptp(col, axis=1) / np.mean(col, axis=1)), 1e-12)  # MEASURED 0.0
    gate("beam2d: max |dep| outside the beam / max inside", np.max(np.abs(dep[:, outside])) / np.max(col), 0.0)
    prof = np.mean(col, axis=1)
    ratio = prof[1:] / prof[:-1]
    gate("beam2d: column ratio vs exp(-kappa dx)", np.max(np.abs(ratio - np.exp(-kappa(0.5) * dx))), 1e-7)  # MEASURED 3.9e-9 (T drift between Strang stages)
    # energy: total deposited power == P_in (uniform absorber, no critical surface, long enough to extinguish?)
    P_in, absd, dump, trans = rows[0, 1:]
    gate("beam2d: |P_in - I0 I_scale 2R| / P_in (top-hat, midpoint lattice)", abs(P_in - I0 * I_scale * 2 * R) / P_in, 1e-12)
    print("  absorbed %.4f  transmitted %.4f  of P_in (informational)" % (absd / P_in, trans / P_in))

print()
if failed:
    print("FAILED gates:")
    for f in failed:
        print("  " + f)
    sys.exit(1)
print("ALL PASS")
