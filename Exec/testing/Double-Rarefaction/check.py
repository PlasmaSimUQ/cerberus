import os
import sys

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

# ---------------------------------------------------------------------------
# Regression check for the 'reconstruction_fallback' option (see
# problem_definition.lua for the problem description and references).
#
# The unlimited O6 reconstruction undershoots the alpha = 0 -> 1 tracer step;
# the minmod fallback repairs the offending faces. At 160x16x16 the final
# min(alpha) is ~ -1.2e-3 without the fallback and ~ -7e-5 with it (the
# residual is update-level leakage: positive faces do not strictly bound the
# updated cell average). rho and p must remain strictly positive throughout.
# ---------------------------------------------------------------------------

ALPHA_MIN_TOL = -3e-4

failed = False

files = get_files(".", include=["plt"], exclude=["temp"], get_all=True)
if not files:
    print("check FAILED: no plotfiles found")
    sys.exit(1)

data = ReadBoxLib(files[-1])

checks = [
    ("alpha_0-fluid", ALPHA_MIN_TOL, "tracer undershoot bounded by fallback"),
    ("rho-fluid", 0.0, "density strictly positive"),
    ("p-fluid", 0.0, "pressure strictly positive"),
]

for name, lower, why in checks:
    x, v = data.get(name)
    ok = v.min() > lower
    print(
        "%s %-14s min = % .6e (require > %g : %s)"
        % ("PASS" if ok else "FAIL", name, v.min(), lower, why)
    )
    failed = failed or not ok

# the fallback must actually have fired (verbosity >= 2 in the lua makes
# calc_reconstruction report per-box hit counts)
n_hits = 0
if os.path.isfile("run_log.txt"):
    with open("run_log.txt") as f:
        n_hits = sum("reconstruction fallback applied" in line for line in f)
ok = n_hits > 0
print(
    "%s fallback activity: %d report lines in run_log.txt (require > 0)"
    % ("PASS" if ok else "FAIL", n_hits)
)
failed = failed or not ok

if failed:
    print("check FAILED")
    sys.exit(1)

print("check PASSED")
