#!/usr/bin/env python3
"""Isotope-scaled stand-in tables from the deuterium eref295 table.

Classical isotope scaling (the H -> D scaling of eos_tools/formats/fpeos.py
and reos3.py, applied D -> T and D -> He4): at equal nuclear NUMBER density
and temperature the per-atom thermodynamics is unchanged, so with
s = m_new / m_D

    rho_new = s rho_D        p_new(rho_new, T)  = p_D(rho_D, T)
    e_new   = e_D / s        cv_new = cv_D / s   (per unit mass)
    dpdrho_new = dpdrho_D / s        dedrho_new = dedrho_D / s^2

On the log10 axes of the .eostab this is a pure relabelling: lrho shifts by
+log10(s), the inverse energy axis le by -log10(s), and the T_of_e / T_of_p /
hull / p / dpdT blocks are copied verbatim. No interpolation. The e_shift of
the common gauge scales by 1/s (per unit mass), which is inconsequential for
the dynamics and ~1e-5 of E_alpha per reaction for the burn bookkeeping.

The He4 stand-in is CRUDE (Z = 1 electronic structure, not helium) and exists
only so the burn product has a table in the mixture. Outputs are
SESAME-derived and must not be committed (data/ is gitignored).
"""

import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "../../python_analysis"))
from eos_tools.formats.eostab import read_eostab, write_eostab  # noqa: E402

AMU_D, AMU_T, AMU_HE4 = 2.014101, 3.016049, 4.002602
src = os.path.join(here, "../EOS-Sesame/data/deuterium_5267_s301_eref295.eostab")
out_dir = os.path.join(here, "data")
os.makedirs(out_dir, exist_ok=True)

prov, meta, blocks = read_eostab(src)
print("source: %s  grid %dx%d  blocks %s" % (os.path.basename(src), len(meta["lrho"]), len(meta["lT"]), sorted(blocks)))

for name, amu, Z in [("tritium", AMU_T, 1.0), ("helium4", AMU_HE4, 2.0)]:
    s = amu / AMU_D
    ls = np.log10(s)
    b = {k: np.array(v, copy=True) for k, v in blocks.items()}
    b["e"] /= s
    b["cv"] /= s
    b["dpdrho"] /= s
    b["dedrho"] /= s * s
    lrho = meta["lrho"] + ls
    inv = {"le": meta["le"] - ls, "lp": np.array(meta["lp"], copy=True)}
    p = [(k, v) for k, v in prov.items() if k not in ("grid", "lrho", "lT", "le", "lp")]
    p = [("material", name + "-standin")] + [(k, v) for k, v in p if k != "material"]
    p = [(k, ("isotope-scaled stand-in from " + v) if k == "source" else v) for k, v in p]
    p = [(k, ("A=%.6f Z=%.1f (isotope-scaled; electronic structure of D) " % (amu, Z)) + v if k == "composition" else v) for k, v in p]
    p = [(k, "%.10e" % (float(v) / s) if k == "e_shift" else v) for k, v in p]
    p.append(("isotope_scaling", "s = m_new/m_D = %.9f: lrho += log10(s), le -= log10(s), e,cv,dpdrho /= s, dedrho /= s^2" % s))
    out = os.path.join(out_dir, "%s_standin_from_D5267_eref295.eostab" % name)
    write_eostab(out, p, lrho, meta["lT"], b, inv)
    # QA: identical per-atom state -> p equal, e per atom equal
    r2, m2, b2 = read_eostab(out)
    i, j = len(lrho) // 2, len(meta["lT"]) // 2
    rho_new, rho_old = 10 ** m2["lrho"][i], 10 ** meta["lrho"][i]
    assert abs(rho_new / rho_old - s) < 1e-9  # axes are written with 12 significant digits
    assert b2["p"][i, j] == blocks["p"][i, j]
    assert abs(b2["e"][i, j] * rho_new / (blocks["e"][i, j] * rho_old) - 1.0) < 1e-9  # energy per volume at equal n
    assert m2["lT"][j] == meta["lT"][j]
    print("  %s: s=%.6f  rho hull [%.3e, %.3e] g/cc  QA ok" % (os.path.basename(out), s, 10 ** m2["lrho"][0], 10 ** m2["lrho"][-1]))
