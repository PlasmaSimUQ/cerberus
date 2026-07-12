#!/usr/bin/env python3
"""Gates for the Dalton-mixture Sod runs (Stage 9 W28 —
doc/eos_mixture_dalton_plan.md section 11).

M1  pure round-off:  mix_pure1 (alpha=1) and mix_pure0 (alpha=0) on two
                     IDENTICAL tables vs the single-table baseline, Linf per
                     field. Every cell takes the pure-cell short-circuit,
                     which is required to be formula-identical to the
                     TabulatedEOS path — any structured difference is a
                     mixture-plumbing bug isolated from the closure math.
M1b uniform mix:     mix_uniform (alpha=0.3) vs single, L1 per field. Dalton
                     is analytically exact for identical ideal-gas
                     components, so the residual measures bilinear
                     interpolation error at the partial densities (the
                     partial-density grid points differ from the total's).
M2  self-test:       MIXEOS-SELFTEST OVERALL PASS in the mix_pure1 (identical
                     tables) and mix_binary (heterogeneous tables) logs.
M3  binary tube:     stable two-material run; alpha in [0,1]; a mixed
                     contact layer of sane width; positivity.
M4  conservation:    total mass, total energy, and PER-COMPONENT mass
                     (sum of the rho*alpha conserved tracer) drift.
M5  wall-time:       mixture advance-time budget vs the single-table run
                     (<= 2.5x for N=2, plan M5). Noise-sensitive at these
                     runtimes — measured values are recorded in the README.
W24 default flux:    no 'flux' key anywhere: the log must show the config
                     default selecting HLLC_general_eos for the mixture gas.

Tolerances marked MEASURED were committed from the first passing run per the
plan's no-guessed-gates rule.
"""

import re
import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

PURE_TOL = 1.0e-11         # Linf/range, formula-identical requirement
UNIFORM_TOL = 1.0e-3       # MEASURED: interpolation-level (same scale as the
                           # EOS-Sod-Ideal twin gate)
CONS_TOL = 1.0e-11         # relative drift
ALPHA_EPS = 1.0e-12        # simplex bound slack
LAYER_MIN, LAYER_MAX = 1, 60   # mixed-cell count bounds (1024 cells, minmod)
WALL_RATIO_MAX = 2.5       # plan M5 budget for N=2

failed = []


def gate(name, ok, detail):
    print("check.py: %-26s %s  %s" % (name, "PASS" if ok else "** FAIL", detail))
    if not ok:
        failed.append(name)


def profile(prefix, which, fields):
    files = sorted(get_files(".", include=[prefix], exclude=["temp"], get_all=True))
    if not files:
        print("check.py: no plotfiles for %s" % prefix)
        sys.exit(1)
    d = ReadBoxLib(files[0 if which == "first" else -1])
    out = {}
    for f in fields:
        x, v = d.get(f)
        if isinstance(x, (list, tuple)):
            x = x[0]
        out[f] = (np.asarray(x).ravel(), np.asarray(v).ravel())
    out["time"] = d.time
    return out


BASE = ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid", "nrg-fluid")
MIX = BASE + ("alpha_0-fluid", "tracer_0-fluid")

single_0 = profile("single.plt", "first", BASE)
single_1 = profile("single.plt", "last", BASE)
pure1_1 = profile("mix_pure1.plt", "last", MIX)
pure0_1 = profile("mix_pure0.plt", "last", MIX)
unif_0 = profile("mix_uniform.plt", "first", MIX)
unif_1 = profile("mix_uniform.plt", "last", MIX)
bin_0 = profile("mix_binary.plt", "first", MIX)
bin_1 = profile("mix_binary.plt", "last", MIX)

# ---------------------------------------------------------------------------
# M1: pure-cell round-off (identical tables, alpha uniform 1 / 0)

for label, d in (("pure1", pure1_1), ("pure0", pure0_1)):
    for f in BASE:
        vs = single_1[f][1]
        vm = d[f][1]
        rng = vs.max() - vs.min()
        linf = np.max(np.abs(vm - vs)) / rng
        gate("%s-roundoff-%s" % (label, f.split("-")[0]), linf <= PURE_TOL,
             "Linf/range=%.3e (tol %g)" % (linf, PURE_TOL))

# ---------------------------------------------------------------------------
# M1b: uniform mixed alpha on identical tables (interpolation-level)

for f in ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid"):
    vs = single_1[f][1]
    vm = unif_1[f][1]
    rng = vs.max() - vs.min()
    l1 = np.mean(np.abs(vm - vs)) / rng
    gate("uniform-twin-%s" % f.split("-")[0], l1 <= UNIFORM_TOL,
         "L1/range=%.3e (tol %g)" % (l1, UNIFORM_TOL))

# ---------------------------------------------------------------------------
# M2: constructor self-test sweeps (grepped from the run logs)

for log in ("run_log_mix_pure1.txt", "run_log_mix_binary.txt"):
    txt = open(log).read()
    gate("selftest-%s" % log.split("_")[-1].split(".")[0],
         "MIXEOS-SELFTEST OVERALL PASS" in txt, "in %s" % log)

# ---------------------------------------------------------------------------
# M3: binary tube — alpha bounds, contact layer, positivity

a = bin_1["alpha_0-fluid"][1]
gate("binary-alpha-bounds", (a.min() >= -ALPHA_EPS) and (a.max() <= 1.0 + ALPHA_EPS),
     "alpha in [%.3e, %.3e]" % (a.min(), a.max()))

n_mixed = int(np.sum((a > 0.01) & (a < 0.99)))
gate("binary-contact-layer", LAYER_MIN <= n_mixed <= LAYER_MAX,
     "%d mixed cells (bounds [%d, %d])" % (n_mixed, LAYER_MIN, LAYER_MAX))

rho_b = bin_1["rho-fluid"][1]
p_b = bin_1["p-fluid"][1]
gate("binary-positivity", (rho_b.min() > 0.0) and (p_b.min() > 0.0),
     "rho_min=%.3e p_min=%.3e" % (rho_b.min(), p_b.min()))

# ---------------------------------------------------------------------------
# M4: conservation — total mass/energy on every run, per-component mass on
#     the binary run (uniform grid: sums are integrals up to a factor)

for name, d0, d1 in (("single", single_0, single_1), ("uniform", unif_0, unif_1),
                     ("binary", bin_0, bin_1)):
    for f, label in (("rho-fluid", "mass"), ("nrg-fluid", "energy")):
        s0 = d0[f][1].sum()
        s1 = d1[f][1].sum()
        drift = abs(s1 - s0) / abs(s0)
        gate("conserve-%s-%s" % (name, label), drift <= CONS_TOL, "drift=%.3e" % drift)

# component 0: the rho*alpha conserved tracer; component 1: rho - tracer
t0 = bin_0["tracer_0-fluid"][1].sum()
t1 = bin_1["tracer_0-fluid"][1].sum()
gate("conserve-binary-comp0", abs(t1 - t0) / abs(t0) <= CONS_TOL,
     "drift=%.3e" % (abs(t1 - t0) / abs(t0)))
c0 = bin_0["rho-fluid"][1].sum() - t0
c1 = bin_1["rho-fluid"][1].sum() - t1
gate("conserve-binary-comp1", abs(c1 - c0) / abs(c0) <= CONS_TOL,
     "drift=%.3e" % (abs(c1 - c0) / abs(c0)))

# ---------------------------------------------------------------------------
# M5: wall-time budget (advance time isolates stepping cost from table parse)


def wall(log):
    txt = open(log).read()
    m = re.findall(r"Run Time advance\s*=\s*([0-9.eE+-]+)", txt)
    return float(m[-1]) if m else None


w_s = wall("run_log_single.txt")
w_b = wall("run_log_mix_binary.txt")
if w_s and w_b:
    ratio = w_b / w_s
    gate("wall-mixture", ratio <= WALL_RATIO_MAX,
         "binary/single = %.2f (%.2fs / %.2fs, max %g)" % (ratio, w_b, w_s, WALL_RATIO_MAX))
else:
    gate("wall-mixture", False, "could not parse 'Run Time advance' from logs")

# ---------------------------------------------------------------------------
# W24 default flux: no run carries a 'flux' key, so the config default must
# have selected the general-EOS solver for the mixture gas

log_txt = open("run_log_mix_binary.txt").read()
gate("default-flux", "defaulting to 'HLLC_general_eos'" in log_txt,
     "config-default notice in run_log_mix_binary.txt")

print("check.py:", "FAIL" if failed else "PASS")
sys.exit(1 if failed else 0)
