#!/usr/bin/env python3
"""Compare the no_fusion / fusion / fusion_brems runs: peak temperatures and
density, total energy and cumulative burn vs time (level-0 data, quarter
domain). Writes compare.png and prints a table."""

import os
import sys

import numpy as np

cmd_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../vis")
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

if not hasattr(np, "product"):
    np.product = np.prod

import matplotlib

matplotlib.use("Agg")
import pylab as plt

m_He4 = 1.987289
CASES = ["no_fusion", "fusion", "fusion_brems"]
here = os.path.dirname(os.path.abspath(__file__))


def history(case):
    d = os.path.join(here, case)
    files = get_files(d, include=[case + ".plt"], exclude=[".inputs", ".txt"], get_all=True)
    out = {k: [] for k in ["t", "Ti", "Te", "rho", "E", "burn", "nD"]}
    nD0 = None
    for f in files:
        r = ReadBoxLib(f, max_level=0)
        (x, y), rho = r.get("rho-ion")
        _, aD = r.get("alpha_0-ion")
        _, aT = r.get("alpha_1-ion")
        _, Ti = r.get("T-ion")
        _, Te = r.get("T-electron")
        _, Ei = r.get("nrg-ion")
        _, Ee = r.get("nrg-electron")
        dV = (x[1] - x[0]) * (y[1] - y[0])
        he = np.clip(1.0 - aD - aT, 0.0, None)
        N = np.sum(rho * he) * dV / m_He4
        nD = np.sum(rho * aD) * dV
        if nD0 is None:
            nD0 = nD
        out["t"].append(r.time)
        out["Ti"].append(Ti.max())
        out["Te"].append(Te.max())
        out["rho"].append(rho.max())
        out["E"].append(np.sum(Ei + Ee) * dV)
        out["burn"].append(N / nD0)
        out["nD"].append(nD)
    return {k: np.array(v) for k, v in out.items()}


H = {c: history(c) for c in CASES}

print("%-13s %8s %8s %8s %10s %10s" % ("case", "max T_i", "max T_e", "max rho", "E_end/E_0", "burn frac"))
for c in CASES:
    h = H[c]
    print("%-13s %8.2f %8.2f %8.2f %10.4f %10.3e" % (c, h["Ti"].max(), h["Te"].max(), h["rho"].max(), h["E"][-1] / h["E"][0], h["burn"][-1]))

fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for c in CASES:
    h = H[c]
    axes[0, 0].plot(h["t"], h["Ti"], label=c)
    axes[0, 1].plot(h["t"], h["Te"], label=c)
    axes[1, 0].plot(h["t"], h["E"] / h["E"][0], label=c)
    axes[1, 1].plot(h["t"], h["burn"], label=c)
for ax, lab in zip(axes.ravel(), [r"max $T_i$ (code)", r"max $T_e$ (code)", r"total fluid energy / $E_0$", "burn fraction"]):
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(lab)
    ax.grid(alpha=0.3)
axes[0, 0].legend()
fig.tight_layout()
fig.savefig(os.path.join(here, "compare.png"), dpi=150)
print("wrote compare.png")
