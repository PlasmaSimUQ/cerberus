#!/usr/bin/env bash

# Common-energy-reference table set (eref295): solid Ti + interior D2 +
# air, one shared energy gauge.
#
#   step (i)  per material: e -> e - e(rho_fill, 295 K)   [--e-ref-state]
#   step (ii) ONE shared positivity constant S            [--e-shift]
#
# S = max over the set of the post-reference full-array deficit + MARGIN,
# so every hull minimum lands O(1)-O(10) code units above zero. Do NOT
# minimise S (several solver paths assume e=0 sits comfortably below
# every hull) and do not inflate it (Newton tolerance is relative to
# |e|).
#
# Run from Exec/testing/EOS-Sesame with the cerberus_python env active and
# the SESAME library available:
#   SES=/path/to/sesame-unc.ascii2 bash make_eref295_set.sh
#
# "One batch, one generator revision": both passes run from the same
# working tree in one invocation; the generator git sha is recorded in
# every header. Outputs are NOT committed (SESAME-derived).

set -euo pipefail

SES=${SES:-/mnt/c/Users/ktap0992/Downloads/sesame-unc/sesame-unc.ascii2}
PREP="python3 ../../python_analysis/eos_table_prep.py"
MARGIN=2.0e13   # erg/g ~ 4.8 code units (e_ref = 4.1293e12 erg/g)

# Reference states: each material's fill density at the 295 K operating
# point; air/D2 get the MOLECULAR floor mass.
TI_ARGS="--mat 2963 --lT 1.25 9.0 384 --lrho -5.30 1.69897 576 \
  --c-cav 5.34 --cold-extend --e-ref-state 4.1856 295"
D2_ARGS="--mat 5267 --table 301 --n-rho 576 --n-T 384 \
  --floor-mass-amu 4.028 --e-ref-state 0.14775 295"
AIR_ARGS="--mat 5031 --table 301 --n-rho 576 --n-T 384 \
  --floor-mass-amu 28.97 --e-ref-state 2.361e-4 295"

TI_OUT=data/ti-beta-21s_2963_coldext_eref295.eostab
D2_OUT=data/deuterium_5267_s301_eref295.eostab
AIR_OUT=data/dry-air_5031_s301_eref295.eostab

probe () {  # run a probe pass, echo the full-array post-reference minimum
    $PREP sesame --src "$SES" $1 --probe-shift \
        | tee /dev/stderr | sed -n 's/.*PROBE-SHIFT min_full=\([^ ]*\).*/\1/p'
}

echo "== pass 1/2: probe post-reference minima =="
MIN_TI=$(probe "$TI_ARGS")
MIN_D2=$(probe "$D2_ARGS")
MIN_AIR=$(probe "$AIR_ARGS")
echo "post-ref minima (erg/g): Ti=$MIN_TI D2=$MIN_D2 air=$MIN_AIR"

S=$(python3 -c "
mins = [float(x) for x in ('$MIN_TI', '$MIN_D2', '$MIN_AIR')]
print('%.10e' % (max(0.0, *(-m for m in mins)) + $MARGIN))")
echo "shared shift S = $S erg/g ($(python3 -c "print('%.3f' % (float('$S')/4.1293e12))") code)"

echo "== pass 2/2: emit the set with the shared gauge =="
$PREP sesame --src "$SES" $TI_ARGS  --e-shift "$S" --out $TI_OUT  --qa qa
$PREP sesame --src "$SES" $D2_ARGS  --e-shift "$S" --out $D2_OUT  --qa qa
$PREP sesame --src "$SES" $AIR_ARGS --e-shift "$S" --out $AIR_OUT --qa qa

echo "== acceptance: design-point gauge spread (target < 0.15 code) =="
python3 - "$S" $TI_OUT $D2_OUT $AIR_OUT <<'EOF'
import sys
sys.path.insert(0, "../../python_analysis")
from eos_tools.condition import sample_bilinear
from eos_tools.formats.eostab import read_eostab

S = float(sys.argv[1])
E_REF_CODE = 4.1293e12  # erg/g per code unit (ref_temp=1e5, ref_mass=m_D)
refs = {"ti": 4.1856, "d2": 0.14775, "air": 2.361e-4}
vals = {}
for tag, path in zip(("ti", "d2", "air"), sys.argv[2:5]):
    prov, meta, blocks = read_eostab(path)
    vals[tag] = sample_bilinear(meta["lrho"], meta["lT"], blocks["e"],
                                refs[tag], 295.0)
    print("%-4s e(design point) = %.10e erg/g  (S deviation %+.3e)"
          % (tag, vals[tag], vals[tag] - S))
spread = (max(vals.values()) - min(vals.values())) / E_REF_CODE
print("design-point gauge spread = %.6e code units (target < 0.15): %s"
      % (spread, "PASS" if spread < 0.15 else "FAIL"))
sys.exit(0 if spread < 0.15 else 1)
EOF

gzip -kf $TI_OUT $D2_OUT $AIR_OUT
echo "done: $TI_OUT $D2_OUT $AIR_OUT (+.gz)"
