#!/usr/bin/env bash

# eref295 set extension: solid aluminum (3720, Ti-style cold extension)
# and solid diamond (7834, native-hull Track-S raise), emitted on the SAME
# shared energy gauge as the original solid/interior/air set
# (make_eref295_set.sh) so any mixture subset of the five tables shares
# one zero of energy.
#
# The gauge is FIXED here, not recomputed: S was chosen by the original
# set's probe pass and is already baked into the shipped D2/air tables.
# New members must adopt it verbatim; the only freedom left is the gate
# that their post-reference minima still land in the O(1)-O(10) code-unit
# band above zero (they do: Al 4.840, diamond 4.846).
#
# Aluminum construction notes (measured, 2026-08-11):
#   - Vinet fit window --fit-rho 2.2 8.0 g/cc: the default (4, 12) is the
#     2963-measured window and extrapolates at Al's rho0=2.70 (ambient
#     residual +3.99 GPa, gate FAIL). The lo edge must include the 306's
#     tension foot (~0.81*rho0) or v0 is unpinned; (2.2, 8.0) minimises
#     rms (1.06e-1) with ambient residual -0.31 GPa.
#   - --z-cond 3: aluminum valence for the Sommerfeld term.
# Diamond has no 411 melt table, so no cold extension is possible; it
# gets the D2/air treatment (native hull, resolution raise). Its T floor
# (72.5 K) narrows diamond-retaining mixture brackets, same class as
# air's 100 K rail.
#
# c_cav is PINNED for both (Ti convention: the T-floor isotherm slope at
# the first density column above the crossover tension foot, ~1.1-1.2x
# rho0 — 2963's 5.34 km/s at 5.8 g/cc). The data-derived default measures
# ON the repaired foot and lands ~2x (Al) / ~11x (diamond) low:
#   Al      5.82 km/s at 3.04 g/cc   (default gave 2.697)
#   diamond 12.7 km/s at 3.83 g/cc   (default gave 1.055)
#
# Run from Exec/testing/EOS-Sesame with the cerberus_python env active:
#   SES=/path/to/sesame-unc.ascii2 bash make_eref295_solids.sh

set -euo pipefail

SES=${SES:-/mnt/c/Users/ktap0992/Downloads/sesame-unc/sesame-unc.ascii2}
PREP="python3 ../../python_analysis/eos_table_prep.py"
S=2.0017551224e13   # erg/g — the original set's shared shift, verbatim
E_REF_CODE=4.1293e12

# Integrity check: if the original set's D2 table is present, its header
# must carry the same S (guards against emitting onto a stale gauge).
D2_REF=data/deuterium_5267_s301_eref295.eostab
if [ -f "$D2_REF" ]; then
    grep -q "e_shift: 2.0017551224e+13" "$D2_REF" \
        || { echo "FATAL: $D2_REF carries a different e_shift than S=$S"; exit 1; }
fi

FOOT_BAR=${FOOT_BAR:-1.0}  # quiet-foot target pressure (bar)

# lrho hi is chosen so log10(rho0 = 2.70) is EXACTLY a lattice node
# (node 493 of 576), as for Ti in make_eref295_set.sh: the reader then
# returns the quiet foot at rho0 itself. rho_hi = 26.79 g/cc (~9.9x
# rho0, was 27).
AL_ARGS="--mat 3720 --lT 1.25 9.0 384 --lrho -5.56 1.4278989136 576 \
  --cold-extend --fit-rho 2.2 8.0 --z-cond 3 --c-cav 5.82 \
  --p-foot $FOOT_BAR --e-ref-state 2.70 295"
C_ARGS="--mat 7834 --table 301 --n-rho 576 --n-T 384 --c-cav 12.7 \
  --p-foot $FOOT_BAR --e-ref-state 3.515 295"

AL_OUT=data/aluminum_3720_coldext_eref295.eostab
C_OUT=data/diamond_7834_s301_eref295.eostab

probe () {  # probe pass, echo the full-array post-reference minimum
    $PREP sesame --src "$SES" $1 --probe-shift \
        | tee /dev/stderr | sed -n 's/.*PROBE-SHIFT min_full=\([^ ]*\).*/\1/p'
}

echo "== pass 1/2: probe post-reference minima against the fixed S =="
MIN_AL=$(probe "$AL_ARGS")
MIN_C=$(probe "$C_ARGS")
python3 - "$S" "$MIN_AL" "$MIN_C" <<'EOF'
import sys
S, e_ref = float(sys.argv[1]), 4.1293e12
ok = True
for tag, m in zip(("Al", "diamond"), sys.argv[2:4]):
    code = (S + float(m)) / e_ref
    band = 0.5 <= code <= 10.0
    ok &= band
    print("%-8s min_final = %.3f code units (band [0.5, 10]): %s"
          % (tag, code, "PASS" if band else "FAIL"))
sys.exit(0 if ok else 1)
EOF

echo "== pass 2/2: emit on the shared gauge =="
$PREP sesame --src "$SES" $AL_ARGS --e-shift "$S" --out $AL_OUT --qa qa
$PREP sesame --src "$SES" $C_ARGS  --e-shift "$S" --out $C_OUT  --qa qa

echo "== acceptance: design-point gauge spread across the full set =="
python3 - "$S" <<'EOF'
import os, sys
sys.path.insert(0, "../../python_analysis")
from eos_tools.condition import sample_bilinear
from eos_tools.formats.eostab import read_eostab

S = float(sys.argv[1])
E_REF_CODE = 4.1293e12
members = [  # (tag, path, ref rho g/cc) — ref T is 295 K for all
    ("ti",  "data/ti-beta-21s_2963_coldext_eref295.eostab", 4.1856),
    ("d2",  "data/deuterium_5267_s301_eref295.eostab", 0.14775),
    ("air", "data/dry-air_5031_s301_eref295.eostab", 2.361e-4),
    ("al",  "data/aluminum_3720_coldext_eref295.eostab", 2.70),
    ("c",   "data/diamond_7834_s301_eref295.eostab", 3.515),
]
vals = {}
for tag, path, rho_ref in members:
    if not os.path.exists(path):
        print("%-4s SKIPPED (not present): %s" % (tag, path)); continue
    prov, meta, blocks = read_eostab(path)
    vals[tag] = sample_bilinear(meta["lrho"], meta["lT"], blocks["e"],
                                rho_ref, 295.0)
    print("%-4s e(design point) = %.10e erg/g  (S deviation %+.3e)"
          % (tag, vals[tag], vals[tag] - S))
spread = (max(vals.values()) - min(vals.values())) / E_REF_CODE
print("design-point gauge spread = %.6e code units (target < 0.15): %s"
      % (spread, "PASS" if spread < 0.15 else "FAIL"))
sys.exit(0 if spread < 0.15 else 1)
EOF

gzip -kf $AL_OUT $C_OUT
echo "done: $AL_OUT $C_OUT (+.gz)"
