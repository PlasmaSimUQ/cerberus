#!/usr/bin/env python3
"""G-e physics check (REPORT-ONLY): the sub-floor rows of the 5031 table vs
the Lemmon et al. 2000 real-gas EOS for this material, as implemented by
CoolProp's pseudo-pure fluid "Air" (validity 59.75-2000 K, <= 2000 MPa; see
doc/eos_air_lowT_extension_plan.md).

For each isotherm in T_CHECK and each table density column up to RHO_MAX,
compare the emitted table against CoolProp:
  * pressure directly;
  * energy as the DROP from the table's own 100 K row, e(100 K) - e(T),
    against CoolProp's u(100 K) - u(T) at the same density (the two energy
    zeros differ; the drop is gauge-free).
CoolProp states that fail (inside the two-phase dome the pseudo-pure
two-phase evaluation may refuse, or the state is off-range) are counted and
skipped. Output: per-isotherm summary (vapor-side cells agree to ~1e-3
by construction — both are the ideal gas; liquid-side cells show the
construction's honest mismatch) plus a CSV in qa/ for plotting.

Usage (from this directory, cerberus_python env):
  python3 qa_coolprop_check.py [data/dry-air_5031_s301_eref295_Tf1p25.eostab]
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../python_analysis"))
from eos_tools.condition import sample_bilinear  # noqa: E402
from eos_tools.formats.eostab import read_eostab  # noqa: E402

T_CHECK = (60.0, 70.0, 80.0, 90.0, 100.0)
RHO_MAX = 0.5          # g/cc — Lemmon's 2 GPa ceiling is ~1 g/cc at these T
VAPOR_SPLIT = 1e-3     # g/cc — report vapor side and dense side separately


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "data/dry-air_5031_s301_eref295_Tf1p25.eostab"
    try:
        import CoolProp.CoolProp as CP
    except ImportError:
        print("CoolProp not importable — G-e report skipped (not a gate)")
        return 0
    prov, meta, blocks = read_eostab(path)
    lrho, lT = meta["lrho"], meta["lT"]
    rho = 10.0 ** lrho
    cols = rho <= RHO_MAX
    rows = []
    print("%-6s %-6s %6s %6s  %-22s %-22s  %s"
          % ("T[K]", "side", "n_ok", "n_fail", "p rel err (median/max)",
             "de rel err (median/max)", "note"))
    for T in T_CHECK:
        stats = {"vapor": [], "dense": []}
        nfail = {"vapor": 0, "dense": 0}
        for r in rho[cols]:
            side = "vapor" if r < VAPOR_SPLIT else "dense"
            try:
                p_cp = CP.PropsSI("P", "Dmass", r * 1e3, "T", T, "Air") * 10.0  # Pa->cgs
                u_cp = CP.PropsSI("Umass", "Dmass", r * 1e3, "T", T, "Air") * 1e4
                u_cp0 = CP.PropsSI("Umass", "Dmass", r * 1e3, "T", 100.0, "Air") * 1e4
                ph = CP.PhaseSI("Dmass", r * 1e3, "T", T, "Air")
            except Exception:
                nfail[side] += 1
                continue
            p_tab = sample_bilinear(lrho, lT, blocks["p"], r, T)
            e_tab = sample_bilinear(lrho, lT, blocks["e"], r, T)
            e_tab0 = sample_bilinear(lrho, lT, blocks["e"], r, 100.0)
            de_tab, de_cp = e_tab0 - e_tab, u_cp0 - u_cp
            rp = abs(p_tab - p_cp) / max(abs(p_cp), 1e-300)
            rde = (abs(de_tab - de_cp) / max(abs(de_cp), 1e-300)) if T < 100.0 else 0.0
            stats[side].append((rp, rde))
            rows.append((T, r, ph, p_cp, p_tab, de_cp, de_tab, rp, rde))
        for side in ("vapor", "dense"):
            s = np.array(stats[side]) if stats[side] else np.zeros((0, 2))
            note = ""
            if side == "vapor" and len(s):
                note = "ideal gas both sides" if np.median(s[:, 0]) < 1e-2 else "CHECK"
            print("%-6g %-6s %6d %6d  %-22s %-22s  %s"
                  % (T, side, len(s), nfail[side],
                     ("%.2e / %.2e" % (np.median(s[:, 0]), s[:, 0].max())) if len(s) else "-",
                     ("%.2e / %.2e" % (np.median(s[:, 1]), s[:, 1].max())) if len(s) else "-",
                     note))
    os.makedirs("qa", exist_ok=True)
    out = os.path.join("qa", "coolprop_check_%s.csv"
                       % os.path.splitext(os.path.basename(path))[0])
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["T_K", "rho_gcc", "coolprop_phase", "p_coolprop_cgs",
                    "p_table_cgs", "de_coolprop", "de_table", "p_rel_err",
                    "de_rel_err"])
        w.writerows(rows)
    print("wrote %s (%d states)" % (out, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
