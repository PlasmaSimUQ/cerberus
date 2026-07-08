#!/usr/bin/env python3
"""Gates for the FPEOS strong-shock Hugoniot case (Stage 4, W9).

1. pre-shock:     the undisturbed target state in each run must satisfy
                  p = p_table(rho, T) — a unit-chain consistency check.
2. on-Hugoniot:   each measured post-shock plateau (rho2, p2) must lie on
                  the Rankine-Hugoniot locus computed HERE, in Python, from
                  the same .eostab file and the run's own pre-shock state:
                    e2(rho2,T2) - e1 = (p2 + p1)/2 * (1/rho1 - 1/rho2)
                  This isolates solver/EOS-path error: the code and the
                  checker share the table but nothing else.
3. physics anchor: the Stage-1 offline locus (qa/hugoniot_D_fpeos.txt,
                  centred at the published cryogenic initial state) must
                  still match the published PIMC points — chains gate 2 to
                  the literature and guards the committed table artifact.
4. conservation:  mass/energy drift, first vs last plotfile, each run.
5. robustness:    the abusive rarefaction run completes (no abort), its
                  hull-clamp counter fired (>0), and its final state is
                  finite everywhere — the W8 machinery is load-bearing.

Tolerances marked MEASURED were committed from the first passing run per
the plan's no-guessed-gates amendment.
"""

import re
import sys

import numpy as np

cmd_folder = "../../../vis"
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

from get_boxlib import ReadBoxLib, get_files

TABLE = "../EOS-Table/data/D_fpeos.eostab"
QA_LOCUS = "../EOS-Table/qa/hugoniot_D_fpeos.txt"
PIMC = "../EOS-Table/data/raw/hugoniot_MC2000_PRL85_1890.txt"

DRIVERS = [4, 40, 400]  # T_DRIVER values run by `run` (code units, x1e5 K)

# --- reference quantities (must mirror problem_definition.lua), in CGS ---
KB = 1.380649e-16          # erg/K
M_D = 3.34358e-24          # g (deuteron)
RHO_REF = 0.171            # g/cc   (171 kg/m^3)
T_REF = 1.0e5              # K
U2_REF = KB * T_REF / M_D  # cm^2/s^2 (u_ref^2)
P_REF = RHO_REF * U2_REF   # barye
GPA = 1.0e-10              # barye -> GPa

PRESHOCK_TOL = 1.0e-5      # MEASURED 2026-07-08: 3.5e-7 (static state, exact)
HUGONIOT_TOL = 0.02        # MEASURED: 0.71/1.07/1.26% (T4/T40/T400), 1024 cells
ANCHOR_TOL = 0.005         # MEASURED: 0.175%; Stage-1 offline 0.08-0.26%
CONS_TOL = 1.0e-11         # MEASURED: <= 3.7e-12 (all waves in-domain)
OUTFLUX_TOL = 1.0e-3       # MEASURED 2026-07-08: 1.9e-15 (undisturbed boundary)
MIN_PLATEAU_CELLS = 8      # MEASURED: 31-43 cells at the chosen stop times

failed = []


def gate(name, ok, detail):
    print("check.py: %-26s %s  %s" % (name, "PASS" if ok else "** FAIL", detail))
    if not ok:
        failed.append(name)


# ---------------------------------------------------------------------------
# minimal .eostab reader (mirrors the frozen spec in EOS-Table/README.md);
# only the p, e and hull blocks are needed here


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
    """Bilinear on (log10 rho, log10 T), clamped to the grid — the same
    evaluation rule as the C++ EosTable."""
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


# ---------------------------------------------------------------------------
# Rankine-Hugoniot locus from the table, centred at (rho1, p1, e1)


def hugoniot_locus(tab, rho1, T1):
    """Return (rho2, p2) arrays [cgs]: for each table temperature above T1,
    the density where the Hugoniot energy equation balances (compressive
    branch, found by sign change on a dense log-density grid)."""
    p1 = bilin(tab, "p", rho1, T1)
    e1 = bilin(tab, "e", rho1, T1)
    # the e block carries a constant e_shift; it cancels in e2 - e1
    rho_grid = np.logspace(np.log10(rho1 * 1.0001), np.log10(rho1 * 12.0), 2000)
    out_r, out_p = [], []
    for lt in tab["lT"]:
        T2 = 10.0 ** lt
        if T2 <= T1 * 1.5:
            continue
        e2 = bilin(tab, "e", rho_grid, T2)
        p2 = bilin(tab, "p", rho_grid, T2)
        f = (e2 - e1) - 0.5 * (p2 + p1) * (1.0 / rho1 - 1.0 / rho_grid)
        s = np.where(np.diff(np.sign(f)) != 0)[0]
        if len(s) == 0:
            continue
        k = s[0]  # first crossing = the physical compressive root
        w = f[k] / (f[k] - f[k + 1])
        r2 = rho_grid[k] * (rho_grid[k + 1] / rho_grid[k]) ** w
        out_r.append(r2)
        out_p.append(float(bilin(tab, "p", r2, T2)))
    return np.asarray(out_r), np.asarray(out_p)


# ---------------------------------------------------------------------------
# plotfile helpers


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
    """Median (rho, p) over the shocked-target plateau — the disturbed
    (p >> p1) cells whose density is near the disturbed maximum. The
    density contrast separates plateau from driver material across the
    contact (shocked cold target ~4x rho1; expanded hot driver much less),
    and from the few smeared shock-ramp cells (medians ignore them)."""
    rho = prof["rho-fluid"][1]
    p = prof["p-fluid"][1]
    rho1 = np.median(rho[-20:])
    p1 = np.median(p[-20:])
    disturbed = p > 2.0 * p1
    rho_max = rho[disturbed].max()
    sel = disturbed & (rho > 0.7 * rho_max)
    return (rho1, p1, np.median(rho[sel]), np.median(p[sel]), int(sel.sum()))


def main():
    # ---------------------------------------------------------------------------
    # 1 + 2: pre-shock consistency and on-Hugoniot, per driver

    tab = load_eostab(TABLE)

    for T in DRIVERS:
        prof = profile("hug_T%d.plt" % T, "last")
        rho1_c, p1_c, rho2_c, p2_c, ncell = plateau(prof)

        # code -> cgs
        rho1, p1 = rho1_c * RHO_REF, p1_c * P_REF
        rho2, p2 = rho2_c * RHO_REF, p2_c * P_REF
        T1 = np.median(prof["T-fluid"][1][-20:]) * T_REF

        # gate 1: the undisturbed state closes through the table
        p1_tab = float(bilin(tab, "p", rho1, T1))
        err1 = abs(p1 - p1_tab) / p1_tab
        gate("preshock-T%d" % T, err1 <= PRESHOCK_TOL and ncell >= MIN_PLATEAU_CELLS,
             "p_sim vs p_table rel err %.2e (tol %g), plateau %d cells"
             % (err1, PRESHOCK_TOL, ncell))

        # gate 2: measured plateau on the table-predicted locus
        locus_r, locus_p = hugoniot_locus(tab, rho1, T1)
        order = np.argsort(locus_p)
        r2_pred = np.interp(np.log(p2), np.log(locus_p[order]), locus_r[order])
        err2 = abs(rho2 - r2_pred) / r2_pred
        gate("hugoniot-T%d" % T, err2 <= HUGONIOT_TOL,
             "compression %.3f vs predicted %.3f at p2=%.1f GPa (rel err %.2e, tol %g)"
             % (rho2 / rho1, r2_pred / rho1, p2 * GPA, err2, HUGONIOT_TOL))

    # ---------------------------------------------------------------------------
    # 3: physics anchor — offline locus vs published PIMC (regression of the
    # committed table artifact; Stage-1 finding: 0.08-0.26% for T >= 250 kK)

    qa = np.loadtxt(QA_LOCUS)                 # compression, P_GPa
    pimc = np.loadtxt(PIMC)                   # T[K], rho, P_GPa, compression
    hot = pimc[pimc[:, 0] >= 2.5e5]
    order = np.argsort(qa[:, 1])
    worst = 0.0
    for Tk, _, P, comp in hot:
        comp_qa = np.interp(np.log(P), np.log(qa[order, 1]), qa[order, 0])
        worst = max(worst, abs(comp_qa - comp) / comp)
    gate("anchor-PIMC", worst <= ANCHOR_TOL,
         "max compression deviation %.2e over %d points (tol %g)"
         % (worst, len(hot), ANCHOR_TOL))

    # ---------------------------------------------------------------------------
    # 4: conservation (uniform grid: sums are integrals up to a constant).
    # Shock runs: stop times keep both waves inside the domain, so the naive
    # sums must be exact. Rarefaction run: the initial |u|=2 outflow drains the
    # boundaries from t=0 by construction; while the inward-moving waves have
    # not reached them the mass loss is exactly the undisturbed boundary flux
    # 2*rho*|u|*t, so the gate checks the drift AGAINST that prediction.

    for T in DRIVERS:
        d0 = profile("hug_T%d.plt" % T, "first")
        d1 = profile("hug_T%d.plt" % T, "last")
        for f, label in (("rho-fluid", "mass"), ("nrg-fluid", "energy")):
            drift = abs(d1[f][1].sum() - d0[f][1].sum()) / abs(d0[f][1].sum())
            gate("conserve-T%d-%s" % (T, label), drift <= CONS_TOL, "drift=%.3e" % drift)

    d0 = profile("rare.plt", "first")
    d1 = profile("rare.plt", "last")
    x = d0["rho-fluid"][0]
    dx = x[1] - x[0]
    m0 = d0["rho-fluid"][1].sum() * dx
    m1 = d1["rho-fluid"][1].sum() * dx
    u_b = abs(d0["x_vel-fluid"][1][0])
    rho_b = d0["rho-fluid"][1][0]
    expect = 2.0 * rho_b * u_b * d1["time"]          # both boundaries drain
    err = abs((m0 - m1) - expect) / expect
    gate("outflux-rarefaction", err <= OUTFLUX_TOL,
         "mass loss %.5f vs boundary-flux prediction %.5f (rel err %.2e, tol %g)"
         % (m0 - m1, expect, err, OUTFLUX_TOL))

    # ---------------------------------------------------------------------------
    # 5: robustness — the abusive run completed, clamped, and stayed finite

    for T in DRIVERS:
        log = open("run_log_T%d.txt" % T).read()
        gate("completed-T%d" % T, "Run Time total" in log, "no abort")

    log = open("run_log_rarefaction.txt").read()
    gate("completed-rarefaction", "Run Time total" in log, "no abort")

    clamps = sum(int(m) for m in re.findall(r"EOS hull clamp applied to (\d+)", log))
    gate("clamps-fired", clamps > 0, "%d clamped evaluations (must be > 0)" % clamps)

    rare = profile("rare.plt", "last")
    finite = all(np.isfinite(rare[f][1]).all()
                 for f in ("rho-fluid", "p-fluid", "x_vel-fluid", "T-fluid", "nrg-fluid"))
    gate("finite-rarefaction", finite, "all fields finite at t_end")

    print("check.py:", "FAIL" if failed else "PASS")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
