#!/usr/bin/env bash

# Sub-floor T extension (doc/eos_air_lowT_extension_plan.md). A mixture is
# bracketed in T by the HIGHEST native floor among its retained components;
# this re-emits the eref295 member carrying that floor (SESAME 5031, native
# floor 100 K) with rows prepended down to lT = 1.25 (~17.7 K, the set's
# shared floor) on the native log spacing, on the SAME shared energy gauge
# as the rest of the eref295 set. Decisions as built:
# D1-a native-spacing prepend, D2-a ideal-thermal branch, D3-a hull 0,
# D4-a --T-floor stage, D5-a NEW file (the shipped table is kept),
# D6-a seam kink accepted.
#
# The gauge is FIXED (S adopted verbatim, cross-checked against another
# member's header) and the extension leaves every native (>= 100 K) node
# untouched, so e(rho_ref, 295 K) and S are unchanged: the other set
# members are NOT re-emitted.
#
# Gates run inline after emission:
#   G-a  native nodes byte-identical to ../data/dry-air_5031_s301_eref295
#        (lrho; lT[k:]; p/e/hull all native rows; dpdrho/dedrho all
#        native rows; dpdT/cv all native rows EXCEPT the seam row, whose
#        PCHIP stencil legitimately changes — reported, not gated)
#   G-b  e_shift header identical; e(ref state) unchanged; 6-member
#        design-point spread < 0.15 code
#   G-c  conditioning record printed (G1 already gated by the CLI)
#
# Run from THIS directory with the cerberus_python env active:
#   SES=/path/to/sesame-unc.ascii2 bash make_air_lowT.sh

set -euo pipefail

SES=${SES:-/mnt/c/Users/ktap0992/Downloads/sesame-unc/sesame-unc.ascii2}
PREP="python3 ../../../python_analysis/eos_table_prep.py"
S=2.0017551224e13   # erg/g — the eref295 shared shift, verbatim
T_FLOOR=${T_FLOOR:-1.25}   # log10 K; the set's shared floor (17.78 K)

REF=../data/dry-air_5031_s301_eref295.eostab
D2_REF=../data/deuterium_5267_s301_eref295.eostab
[ -f "$REF" ] || { echo "FATAL: shipped reference table $REF missing (gunzip it)"; exit 1; }
grep -q "e_shift: 2.0017551224e+13" "$D2_REF" \
    || { echo "FATAL: $D2_REF carries a different e_shift than S=$S"; exit 1; }
grep -q "e_shift: 2.0017551224e+13" "$REF" \
    || { echo "FATAL: $REF carries a different e_shift than S=$S"; exit 1; }

# IDENTICAL to this member's arguments in make_eref295_set.sh + --T-floor
AIR_ARGS="--mat 5031 --table 301 --n-rho 576 --n-T 384 \
  --floor-mass-amu 28.97 --e-ref-state 2.361e-4 295 --T-floor $T_FLOOR"
OUT=data/dry-air_5031_s301_eref295_Tf1p25.eostab

echo "== emit on the shared gauge (S = $S) =="
$PREP sesame --src "$SES" $AIR_ARGS --e-shift "$S" --out $OUT --qa qa

echo "== gates G-a / G-b =="
python3 - "$OUT" "$REF" "$S" <<'PY'
import os, sys
import numpy as np
sys.path.insert(0, "../../../python_analysis")
from eos_tools.condition import sample_bilinear
from eos_tools.formats.eostab import read_eostab

new, ref, S = sys.argv[1], sys.argv[2], float(sys.argv[3])
pn, mn, bn = read_eostab(new)
pr, mr, br = read_eostab(ref)
ok = True
def gate(name, cond):
    global ok
    print("  %-58s %s" % (name, "PASS" if cond else "FAIL"))
    ok &= bool(cond)
k = len(mn["lT"]) - len(mr["lT"])
print("G-a native byte identity (k = %d prepended rows)" % k)
gate("lrho identical", np.array_equal(mn["lrho"], mr["lrho"]))
dlT = np.abs(mn["lT"][k:] - mr["lT"]).max()
gate("lT[k:] identical to <= 1e-11 (header stores 12-digit endpoints; "
     "nodes are reconstructed) max %.2e" % dlT, dlT <= 1e-11)
for blk in ("p", "e", "hull", "dpdrho", "dedrho"):
    gate("%s identical on all native rows" % blk,
         np.array_equal(bn[blk][:, k:], br[blk]))
for blk in ("dpdT", "cv"):
    gate("%s identical on native rows except the seam row" % blk,
         np.array_equal(bn[blk][:, k + 1:], br[blk][:, 1:]))
    a, b = bn[blk][:, k], br[blk][:, 0]
    rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
    print("  %-58s max rel %.3e median %.3e (report)" %
          ("%s seam-row (100 K) stencil change" % blk, rel.max(), np.median(rel)))
gate("new rows are hull 0", np.all(bn["hull"][:, :k] == 0.0))
gate("new rows p > 0 and strictly increasing in T (incl. seam)",
     bn["p"][:, :k].min() > 0 and np.all(np.diff(bn["p"][:, :k + 1], axis=1) > 0))
gate("new rows e strictly increasing in T (incl. seam)",
     np.all(np.diff(bn["e"][:, :k + 1], axis=1) > 0))
print("G-b gauge")
hn = dict(pn); hr = dict(pr)
gate("e_shift header identical (%s)" % hn.get("e_shift"),
     hn.get("e_shift") == hr.get("e_shift") and
     abs(float(hn["e_shift"]) - S) <= 1e-6 * S)
e_new = sample_bilinear(mn["lrho"], mn["lT"], bn["e"], 2.361e-4, 295.0)
e_ref = sample_bilinear(mr["lrho"], mr["lT"], br["e"], 2.361e-4, 295.0)
gate("e(2.361e-4 g/cc, 295 K) unchanged (%.10e vs %.10e)" % (e_new, e_ref),
     abs(e_new - e_ref) <= 1e-11 * abs(e_ref))
gate("min(e) + 0 > 0 (forced S still clears the new rows): min %.4e" % bn["e"].min(),
     bn["e"].min() > 0)
E_REF_CODE = 4.1293e12
members = [("ti", "../data/ti-beta-21s_2963_coldext_eref295.eostab", 4.1856),
           ("d2", "../data/deuterium_5267_s301_eref295.eostab", 0.14775),
           ("air", ref, 2.361e-4), ("air-lowT", new, 2.361e-4),
           ("al", "../data/aluminum_3720_coldext_eref295.eostab", 2.70),
           ("c", "../data/diamond_7834_s301_eref295.eostab", 3.515)]
vals = {}
for tag, path, rho_ref in members:
    if not os.path.exists(path):
        print("  %-8s SKIPPED (not present): %s" % (tag, path)); continue
    _, m, b = read_eostab(path)
    vals[tag] = sample_bilinear(m["lrho"], m["lT"], b["e"], rho_ref, 295.0)
spread = (max(vals.values()) - min(vals.values())) / E_REF_CODE
gate("design-point gauge spread %.6e code (< 0.15) over %d members"
     % (spread, len(vals)), spread < 0.15)
print("G-c conditioning record:\n  %s" % hn.get("conditioning", "?"))
print("T axis: %d nodes, lT [%.6g, %.6g] -> T [%.4g, %.4g] K; hull coverage %.2f%%"
      % (len(mn["lT"]), mn["lT"][0], mn["lT"][-1], 10**mn["lT"][0], 10**mn["lT"][-1],
         100.0 * (bn["hull"] > 0.5).mean()))
sys.exit(0 if ok else 1)
PY

gzip -kf $OUT
echo "done: $OUT (+.gz)"
