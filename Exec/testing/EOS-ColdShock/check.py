#!/usr/bin/env python3
"""Gates for the SS4 cold-start shocks on the spliced deuterium table.

1. pre-shock:     the undisturbed 24 K liquid must be EXACTLY the loaded
                  initial state: rho and T closure, internal-energy
                  closure through the table, and |p_sim - p_table|
                  normalized by the SHOCKED-PLATEAU pressure — the cold
                  reference cell's own p is resolution-sensitive at the
                  tens-of-bar level (doc/eos_creation_plan.md §1 trap),
                  so a relative-p gate there would be meaningless.
2. on-Hugoniot:   each measured post-shock plateau (rho2, p2) must lie on
                  the Rankine-Hugoniot locus computed HERE from the same
                  .eostab and the run's own pre-shock state. The code and
                  the checker share the table but nothing else.
3. no-clamp:      the hull-clamp counter must stay ZERO for the physical
                  runs — the whole path from 24 K liquid to plasma lives
                  inside the table hull (the SS4-specific gate).
4. anchor-PIMC:   the table's own cold-start locus (computed here) must
                  match the published MC2000 PIMC points — chains gate 2
                  to the literature without any offline artifact.
5. conservation:  mass/energy drift, first vs last plotfile, per run;
                  boundary-flux-predicted drift for the rarefaction.
6. robustness:    the abusive dome-expansion run completes, its clamps
                  fire (>0), and its final state is finite everywhere.

Tolerances marked MEASURED were committed from the first passing run per
the house no-guessed-gates convention.
"""

import re
import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

TABLE = "../EOS-Table/data/D_spliced.eostab"
PIMC = "../EOS-Table/data/raw/hugoniot_MC2000_PRL85_1890.txt"

RUNS = ["gasgun", "mbar", "plasma"]

# --- reference quantities (must mirror problem_definition.lua), in CGS ---
KB = 1.380649e-16          # erg/K
M_D = 3.34358e-24          # g (deuteron)
RHO_REF = 0.171            # g/cc   (171 kg/m^3)
T_REF = 1.0e5              # K
U2_REF = KB * T_REF / M_D  # cm^2/s^2 (u_ref^2)
P_REF = RHO_REF * U2_REF   # barye
GPA = 1.0e-10              # barye -> GPa

T_COLD = 24.0              # K, the initial target temperature

PRESHOCK_RHOT_TOL = 1.0e-12  # MEASURED 2026-07-12: <=2e-13 (static init)
PRESHOCK_E_TOL = 1.0e-6      # MEASURED: 3.48e-7, identical in all three
                             # runs (systematic nondim unit-chain
                             # round-off, not state drift)
PRESHOCK_P_TOL = 1.0e-4      # |p_sim-p_tab|/p2; MEASURED: <=3e-10
HUGONIOT_TOL = 0.02          # MEASURED: 0.36/0.85/0.60% (gasgun/mbar/plasma)
ANCHOR_MED_TOL = 0.05        # MEASURED: 1.2% median (SS3 offline: 1.19%)
CONS_TOL = 1.0e-11           # MEASURED: <=4e-12 (all waves in-domain)
OUTFLUX_TOL = 1.0e-3
MIN_PLATEAU_CELLS = 8

failed = []


def gate(name, ok, detail):
    print("check.py: %-26s %s  %s" % (name, "PASS" if ok else "** FAIL", detail))
    if not ok:
        failed.append(name)


# ---------------------------------------------------------------------------
# minimal .eostab reader + bilinear rule (deliberately independent of
# eos_tools — the checker shares only the table file with the generator)


def load_eostab(path):
    hdr = {}
    blocks = {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("block:"):
                cur = line.split(":", 1)[1].strip()
                blocks[cur] = []
                continue
            if cur is None:
                if ":" in line:
                    k, v = line.split(":", 1)
                    hdr[k.strip()] = v.strip()
                continue
            blocks[cur].extend(float(t) for t in line.split())

    n_rho, n_T = (int(w.split("=")[1]) for w in hdr["grid"].split())
    lrho = tuple(float(t) for t in hdr["lrho"].split())
    lT = tuple(float(t) for t in hdr["lT"].split())
    tab = {
        "n_rho": n_rho, "n_T": n_T,
        "lrho": np.linspace(lrho[0], lrho[1], n_rho),
        "lT": np.linspace(lT[0], lT[1], n_T),
    }
    for name in ("p", "e", "hull"):
        tab[name] = np.asarray(blocks[name]).reshape(n_rho, n_T)
    return tab


def bilin(tab, name, rho, T):
    lr = np.log10(rho)
    lt = np.log10(T)
    lrho, lT = tab["lrho"], tab["lT"]
    i = np.clip(np.searchsorted(lrho, lr) - 1, 0, tab["n_rho"] - 2)
    j = np.clip(np.searchsorted(lT, lt) - 1, 0, tab["n_T"] - 2)
    fx = np.clip((lr - lrho[i]) / (lrho[i + 1] - lrho[i]), 0.0, 1.0)
    fy = np.clip((lt - lT[j]) / (lT[j + 1] - lT[j]), 0.0, 1.0)
    v = tab[name]
    return ((1 - fx) * (1 - fy) * v[i, j] + fx * (1 - fy) * v[i + 1, j] +
            (1 - fx) * fy * v[i, j + 1] + fx * fy * v[i + 1, j + 1])


def hugoniot_locus(tab, rho1, T1):
    """(rho2, p2) [cgs] per table temperature: compressive-branch root of
    the RH energy equation. Rows below 10*T1 are skipped — they thread the
    two-phase dome's hull-0 fill, and every measured plateau sits far
    above them."""
    p1 = bilin(tab, "p", rho1, T1)
    e1 = bilin(tab, "e", rho1, T1)
    rho_grid = np.logspace(np.log10(rho1 * 1.0001), np.log10(rho1 * 12.0), 2000)
    out_r, out_p = [], []
    for lt in tab["lT"]:
        T2 = 10.0 ** lt
        if T2 <= 10.0 * T1:
            continue
        e2 = bilin(tab, "e", rho_grid, T2)
        p2 = bilin(tab, "p", rho_grid, T2)
        f = (e2 - e1) - 0.5 * (p2 + p1) * (1.0 / rho1 - 1.0 / rho_grid)
        s = np.where(np.diff(np.sign(f)) != 0)[0]
        if len(s) == 0:
            continue
        k = s[0]
        w = f[k] / (f[k] - f[k + 1])
        r2 = rho_grid[k] * (rho_grid[k + 1] / rho_grid[k]) ** w
        out_r.append(r2)
        out_p.append(float(bilin(tab, "p", r2, T2)))
    return np.asarray(out_r), np.asarray(out_p)


def profile(prefix, which):
    files = sorted(get_files(".", include=[prefix], exclude=["temp"], get_all=True))
    if not files:
        print("check.py: no plotfiles for %s" % prefix)
        sys.exit(1)
    d = ReadBoxLib(files[0 if which == "first" else -1])
    out = {}
    for f in ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid", "nrg-fluid"):
        x, v = d.get(f)
        if isinstance(x, (list, tuple)):
            x = x[0]
        out[f] = (np.asarray(x).ravel(), np.asarray(v).ravel())
    out["time"] = d.time
    return out


def plateau(prof):
    """Median (rho, p) over the shocked-target plateau (see EOS-Hugoniot:
    density contrast separates shocked target from expanded driver)."""
    rho = prof["rho-fluid"][1]
    p = prof["p-fluid"][1]
    rho1 = np.median(rho[-20:])
    p1 = np.median(p[-20:])
    disturbed = p > max(2.0 * p1, 1e-3)
    rho_max = rho[disturbed].max()
    sel = disturbed & (rho > 0.7 * rho_max)
    return (rho1, p1, np.median(rho[sel]), np.median(p[sel]), int(sel.sum()))


def main():
    tab = load_eostab(TABLE)

    for tag in RUNS:
        prof = profile("cs_%s.plt" % tag, "last")
        rho1_c, p1_c, rho2_c, p2_c, ncell = plateau(prof)

        rho1, p1 = rho1_c * RHO_REF, p1_c * P_REF
        rho2, p2 = rho2_c * RHO_REF, p2_c * P_REF
        T1 = np.median(prof["T-fluid"][1][-20:]) * T_REF

        # gate 1: the undisturbed cold state is exactly the loaded state
        drho = abs(rho1_c - 1.0)
        dT = abs(T1 / T_COLD - 1.0)
        e1_sim = np.median(prof["nrg-fluid"][1][-20:] /
                           prof["rho-fluid"][1][-20:]) * U2_REF  # at rest
        e1_tab = float(bilin(tab, "e", rho1, T1))
        de = abs(e1_sim / e1_tab - 1.0)
        p1_tab = float(bilin(tab, "p", rho1, T1))
        dp = abs(p1 - p1_tab) / p2      # plateau-normalized (see docstring)
        gate("preshock-%s" % tag,
             drho <= PRESHOCK_RHOT_TOL and dT <= 1e-6
             and de <= PRESHOCK_E_TOL and dp <= PRESHOCK_P_TOL
             and ncell >= MIN_PLATEAU_CELLS,
             "drho=%.1e dT=%.1e de=%.2e dp/p2=%.2e plateau %d cells"
             % (drho, dT, de, dp, ncell))

        # gate 2: measured plateau on the table-predicted locus
        locus_r, locus_p = hugoniot_locus(tab, rho1, T1)
        order = np.argsort(locus_p)
        r2_pred = np.interp(np.log(p2), np.log(locus_p[order]), locus_r[order])
        err2 = abs(rho2 - r2_pred) / r2_pred
        gate("hugoniot-%s" % tag, err2 <= HUGONIOT_TOL,
             "compression %.3f vs predicted %.3f at p2=%.4g GPa (rel err "
             "%.2e, tol %g)" % (rho2 / rho1, r2_pred / rho1, p2 * GPA,
                                err2, HUGONIOT_TOL))

        # gate 3 (SS4): zero hull clamps on the physical path
        log = open("run_log_%s.txt" % tag).read()
        clamps = sum(int(m) for m in
                     re.findall(r"EOS hull clamp applied to (\d+)", log))
        gate("noclamp-%s" % tag, clamps == 0,
             "%d clamped evaluations (must be 0)" % clamps)
        gate("completed-%s" % tag, "Run Time total" in log, "no abort")

    # gate 4: table's own cold locus vs published PIMC (median deviation;
    # the low-P PIMC points carry their own large error bars)
    pimc = np.loadtxt(PIMC)  # T rho P_GPa compression
    locus_r, locus_p = hugoniot_locus(tab, 0.171, T_COLD)
    order = np.argsort(locus_p)
    devs = []
    for Tk, _, P_gpa, comp in pimc:
        P_b = P_gpa * 1.0e10  # GPa -> barye
        if P_b < locus_p.min() or P_b > locus_p.max():
            continue
        c_tab = np.interp(np.log(P_b), np.log(locus_p[order]),
                          locus_r[order]) / 0.171
        devs.append(abs(c_tab / comp - 1.0))
    devs = np.asarray(devs)
    gate("anchor-PIMC", len(devs) > 0 and float(np.median(devs)) <= ANCHOR_MED_TOL,
         "median %.2f%% max %.2f%% over %d in-range points (tol median %g)"
         % (100 * np.median(devs), 100 * devs.max(), len(devs),
            ANCHOR_MED_TOL))

    # gate 5: conservation
    for tag in RUNS:
        d0 = profile("cs_%s.plt" % tag, "first")
        d1 = profile("cs_%s.plt" % tag, "last")
        for f, label in (("rho-fluid", "mass"), ("nrg-fluid", "energy")):
            drift = abs(d1[f][1].sum() - d0[f][1].sum()) / abs(d0[f][1].sum())
            gate("conserve-%s-%s" % (tag, label), drift <= CONS_TOL,
                 "drift=%.3e" % drift)

    d0 = profile("rare.plt", "first")
    d1 = profile("rare.plt", "last")
    x = d0["rho-fluid"][0]
    dx = x[1] - x[0]
    m0 = d0["rho-fluid"][1].sum() * dx
    m1 = d1["rho-fluid"][1].sum() * dx
    u_b = abs(d0["x_vel-fluid"][1][0])
    rho_b = d0["rho-fluid"][1][0]
    expect = 2.0 * rho_b * u_b * d1["time"]
    err = abs((m0 - m1) - expect) / expect
    gate("outflux-rarefaction", err <= OUTFLUX_TOL,
         "mass loss %.5f vs boundary-flux prediction %.5f (rel err %.2e)"
         % (m0 - m1, expect, err))

    # gate 6: robustness of the abusive dome expansion
    log = open("run_log_rarefaction.txt").read()
    gate("completed-rarefaction", "Run Time total" in log, "no abort")
    clamps = sum(int(m) for m in
                 re.findall(r"EOS hull clamp applied to (\d+)", log))
    gate("clamps-fired", clamps > 0,
         "%d clamped evaluations (must be > 0)" % clamps)
    rare = profile("rare.plt", "last")
    finite = all(np.isfinite(rare[f][1]).all()
                 for f in ("rho-fluid", "p-fluid", "x_vel-fluid",
                           "T-fluid", "nrg-fluid"))
    gate("finite-rarefaction", finite, "all fields finite at t_end")

    print("check.py:", "FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
