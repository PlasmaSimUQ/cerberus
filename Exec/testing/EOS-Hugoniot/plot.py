#!/usr/bin/env python3
"""Hugoniot overlay: simulated post-shock states vs the table-predicted
locus (left), and the offline table locus vs published PIMC (right)."""


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
from check import (DRIVERS, GPA, P_REF, RHO_REF, T_REF, bilin, hugoniot_locus,
                   load_eostab, plateau, profile, TABLE, QA_LOCUS, PIMC)

tab = load_eostab(TABLE)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

# --- left: simulation vs sim-centred locus -------------------------------
prof = profile("hug_T%d.plt" % DRIVERS[0], "last")
rho1_c, _, _, _, _ = plateau(prof)
rho1 = rho1_c * RHO_REF
T1 = np.median(prof["T-fluid"][1][-20:]) * T_REF

lr, lp = hugoniot_locus(tab, rho1, T1)
ax1.plot(lp * GPA, lr / rho1, "-", color="0.35", lw=1.5,
         label="RH locus from table\n(centred at sim initial state)")

for T, mk in zip(DRIVERS, ("o", "s", "^")):
    prof = profile("hug_T%d.plt" % T, "last")
    r1, p1, r2, p2, _ = plateau(prof)
    ax1.plot(p2 * P_REF * GPA, r2 / r1, mk, ms=8, mfc="none", mew=1.8,
             label="simulated, $T_{driver}$=%.0e K" % (T * T_REF))

ax1.set_xscale("log")
ax1.set_xlabel("post-shock pressure [GPa]")
ax1.set_ylabel(r"compression $\rho_2/\rho_1$")
ax1.set_title("Cerberus (effective-$\\gamma$) vs Rankine-Hugoniot + FPEOS table")
ax1.legend(fontsize=8, loc="lower right")
ax1.grid(alpha=0.3)

# --- right: offline anchor (Stage-1 artifact) vs published PIMC ----------
qa = np.loadtxt(QA_LOCUS)
pimc = np.loadtxt(PIMC)
ax2.plot(qa[:, 1], qa[:, 0], "-", color="0.35", lw=1.5,
         label="table locus (offline QA,\ncryogenic initial state)")
ax2.plot(pimc[:, 2], pimc[:, 3], "d", ms=7, mfc="none", mew=1.8, color="C3",
         label="PIMC, Militzer & Ceperley\nPRL 85, 1890 (2000)")
ax2.set_xscale("log")
ax2.set_xlabel("pressure [GPa]")
ax2.set_ylabel(r"compression $\rho/\rho_0$")
ax2.set_title("physics anchor: table vs published locus")
ax2.legend(fontsize=8, loc="lower right")
ax2.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("plot.png", dpi=150)
print("wrote plot.png")
