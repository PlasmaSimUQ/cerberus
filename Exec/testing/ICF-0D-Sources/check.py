#!/usr/bin/env python3
"""Gates for the 0-D ICF source-term checks (see problem_definition.lua).

brems     : with n fixed and E = 1.5 n T, dE/dt = -P_br gives the closed form
            sqrt(T(t)) = sqrt(T0) - a t / 2,  a = C_br n_e Z2n / (1.5 n_ion).
            Also: rho and both fuel tracers unchanged (brems moves no mass).
fusion    : exact mass identities per reaction count N,
              (rho_D0 - rho_D)/m_D = (rho_T0 - rho_T)/m_T = rho_He4/m_He4 = N,
              rho0 - rho = (m_D + m_T - m_He4) N,  nrg unchanged (v = 0, alpha_fraction = 0);
            and n_D(t) against the ODE dN/dt = n_D n_T <sv>(T) with
            T = p0/(n0 - N) (energy fixed, ion count falls), RK4 fine-stepped.
ignition  : T falls at 3.5 keV and rises at 6 keV (P_alpha = P_br at ~4.6 keV
            with the NRL fit and the 1.69e-38 bremsstrahlung coefficient).

References: NRL Plasma Formulary (Huba, NRL, 2019) for the bremsstrahlung
coefficient and the DT reactivity fit; Fraley et al., Phys. Fluids 17, 474
(1974) for the alpha partition (not exercised here: single fluid).

Tolerances marked MEASURED were committed from the first passing run.
"""

import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

if not hasattr(np, "product"):
    np.product = np.prod

# ---- the same constants and reference quantities as the lua / the actions ----
kB, e, c = 1.38064852e-23, 1.602176634e-19, 299792458.0
x_ref, rho_ref, m_ref = 1.0e-4, 1000.0, 3.34358e-27
T_ref = 1.0e3 * e / kB
n_ref = rho_ref / m_ref
u_ref = np.sqrt(kB * T_ref / m_ref)
t_ref = x_ref / u_ref
prs_ref = n_ref * m_ref * u_ref**2
kT_eV = kB * T_ref / e
kT_keV = kT_eV / 1e3

C_br = 1.69e-38 * n_ref**2 * np.sqrt(kT_eV) * t_ref / prs_ref
C_sv = 3.68e-18 * n_ref * t_ref
E_alpha = 3.5e6 * e / (m_ref * u_ref**2)

m_D, m_T, m_He4 = 1.0, 1.497877, 1.987289
m_loss = m_D + m_T - m_He4  # fluid mass lost per reaction (neutron + Q/c^2)
rho0 = 1.0
n0 = rho0 / (0.5 * (m_D + m_T))
alpha_D = m_D / (m_D + m_T)

S = "plasma"
failed = []


def gate(name, value, tol):
    ok = value <= tol
    print("  %-44s %.3e  (tol %.1e)  %s" % (name, value, tol, "PASS" if ok else "FAIL"))
    if not ok:
        failed.append(name)


def series(prefix):
    files = get_files(".", include=[prefix + ".plt"], exclude=[".inputs", ".txt"], get_all=True)
    out = {"t": [], "rho": [], "T": [], "nrg": [], "rho_D": [], "rho_T": []}
    for f in files:
        ds = ReadBoxLib(f, max_level=0)
        out["t"].append(ds.time)
        for key, comp in [("rho", "rho-"), ("T", "T-"), ("nrg", "nrg-"),
                          ("rho_D", "tracer_0-"), ("rho_T", "tracer_1-")]:
            _, d = ds.get(comp + S)
            out[key].append(np.mean(d))  # uniform box: mean == every cell
    return {k: np.array(v) for k, v in out.items()}


def sv(T_keV):
    x = T_keV ** (-1.0 / 3.0)
    return C_sv * x * x * np.exp(-19.94 * x)


# ------------------------------------------------------------------ brems
print("brems")
d = series("brems")
T0 = d["T"][0]
a = C_br * n0 * n0 / (1.5 * n0)  # n_e = Z2n = n_ion for a Z = 1 plasma
T_exact = (np.sqrt(T0) - 0.5 * a * d["t"]) ** 2
gate("brems: max |T - exact| / T0", np.max(np.abs(d["T"] - T_exact)) / T0, 1e-9)  # MEASURED 2.4e-11
gate("brems: max |rho - rho0|", np.max(np.abs(d["rho"] - rho0)), 1e-12)
gate("brems: max |rho_D - rho_D0|", np.max(np.abs(d["rho_D"] - d["rho_D"][0])), 1e-12)
gate("brems: max |rho_T - rho_T0|", np.max(np.abs(d["rho_T"] - d["rho_T"][0])), 1e-12)

# ------------------------------------------------------------------ fusion
print("fusion")
d = series("fusion")
rho_He4 = d["rho"] - d["rho_D"] - d["rho_T"]
N = (d["rho_D"][0] - d["rho_D"]) / m_D  # reactions per volume from D
Nend = N[-1]
gate("fusion: |dT/m_T - N| / N", np.max(np.abs((d["rho_T"][0] - d["rho_T"]) / m_T - N)[1:] / N[1:]), 1e-10)
gate("fusion: |He4/m_He4 - N| / N", np.max(np.abs(rho_He4 / m_He4 - N)[1:] / N[1:]), 1e-10)
gate("fusion: |(rho0 - rho) - m_loss N| / (m_loss N)", np.max(np.abs((rho0 - d["rho"]) - m_loss * N)[1:] / (m_loss * N[1:])), 1e-10)
gate("fusion: max |nrg - nrg0| / nrg0", np.max(np.abs(d["nrg"] - d["nrg"][0])) / d["nrg"][0], 1e-12)

# ODE for N(t): dN/dt = (nD0 - N)(nT0 - N) <sv>(p0 / (n0 - N))
nD0 = nT0 = 0.5 * n0
p0 = n0 * d["T"][0]


def rhs(N):
    return (nD0 - N) * (nT0 - N) * sv(p0 / (n0 - N) * kT_keV)


t_out, N_num = d["t"], 0.0
N_ode, t = [0.0], 0.0
h = (t_out[-1] - t_out[0]) / 200000.0
for tk in t_out[1:]:
    while t < tk - 1e-12:
        hh = min(h, tk - t)
        k1 = rhs(N_num)
        k2 = rhs(N_num + 0.5 * hh * k1)
        k3 = rhs(N_num + 0.5 * hh * k2)
        k4 = rhs(N_num + hh * k3)
        N_num += hh / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += hh
    N_ode.append(N_num)
N_ode = np.array(N_ode)
gate("fusion: max |N - ode| / N_end", np.max(np.abs(N - N_ode)) / Nend, 1e-9)  # MEASURED 2.5e-11
print("  burn fraction over the run: %.3e" % (Nend / nD0))

# ------------------------------------------------------------------ ignition
print("ignition")
for tag, sign in [("ignition_lo", -1.0), ("ignition_hi", +1.0)]:
    d = series(tag)
    dT = d["T"][-1] - d["T"][0]
    print("  %s: T0 = %.2f keV, dT = %+.4e" % (tag, d["T"][0], dT))
    if sign * dT <= 0.0:
        failed.append(tag + ": wrong sign of dT")
        print("  FAIL")
    else:
        print("  PASS")

print()
if failed:
    print("FAILED gates:")
    for f in failed:
        print("  " + f)
    sys.exit(1)
print("ALL PASS")
