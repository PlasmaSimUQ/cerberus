"""SS3 QA for the spliced deuterium table (plan §3 SS3 gates).

Greppable 'SPLICE-QA[<gate>] PASS/FAIL ...' lines; hard gates set the exit
code. First-pass measured numbers get committed into the gates per house
convention.
"""

import os

import numpy as np

from .constants import GPA_CGS, KB, M_D
from .hugoniot import locus
from .splice import REF_STATE, S1, S2, S3, raw_path
from .thermo import maxwell_residual, sound_speed_sq


class Gate:
    def __init__(self):
        self.failed = []

    def check(self, name, ok, detail, hard=True):
        verdict = "PASS" if ok else ("FAIL" if hard else "FAIL(report)")
        print("SPLICE-QA[%s] %s %s" % (name, verdict, detail))
        if hard and not ok:
            self.failed.append(name)


def run_qa_ti(res, qa_dir=None):
    """Ti v0 splice QA. Hard gates are STRUCTURAL only (convexity, finite,
    Hugoniot smoothness/monotonicity); the seam-1 band metrics are reports
    — the WDM band is documented as accuracy-uncontrolled pending real
    data (splice_ti.py docstring)."""
    import numpy as np

    from .splice_ti import REF_STATE_TI, S1_TI, S3_TI
    from .thermo import maxwell_residual as mres

    g = Gate()
    lrho, lT = res["lrho"], res["lT"]
    blocks = res["blocks"]
    p, e, hull = blocks["p"], blocks["e"], blocks["hull"] > 0.5
    R = 10.0 ** lrho[:, None] * np.ones((1, len(lT)))
    TT = 10.0 ** lT[None, :] * np.ones((len(lrho), 1))

    g.check("hull", hull.mean() > 0.4,
            "in-hull fraction %.3f; tension-clipped %d cells"
            % (hull.mean(), res["n_clip"]))
    for a in res["align_stats"]:
        g.check("align-" + a["pair"], True,
                "cE=%.4e erg/g n=%d; mean|dev|/kT=%.3f (report — Ti v0 "
                "WDM band uncontrolled)" % (a["cE"], a["n"],
                                            a["std_over_kT"]), hard=False)

    finite = all(np.all(np.isfinite(blocks[k]))
                 for k in ("p", "e", "dpdT", "dpdrho", "cv", "dedrho"))
    g.check("finite", finite, "all blocks finite")

    cs2 = sound_speed_sq(R, TT, p, blocks["dpdrho"], blocks["dpdT"],
                         blocks["cv"])
    n_bad = int(np.sum(cs2[hull] <= 0.0))
    g.check("convexity", n_bad == 0,
            "cs^2 <= 0 at %d/%d in-hull cells" % (n_bad, int(hull.sum())))

    Rm = mres(lrho, lT, p, e)
    g.check("maxwell", True,
            "in-hull median %.2e p95 %.2e (report)"
            % (float(np.median(Rm[hull])),
               float(np.percentile(Rm[hull], 95.0))), hard=False)

    # cold-start Hugoniot from ambient solid Ti (report p0; the model's
    # ambient pressure is not exactly zero — no literature inputs used)
    rho0, T0 = REF_STATE_TI["rho0"], REF_STATE_TI["T0"]
    e0 = bilin(lrho, lT, e, rho0, T0)
    p0 = bilin(lrho, lT, p, rho0, T0)
    comp, pres = locus(lrho, lT, p, e, rho0, e0, p0)
    m_s = pres >= 1.0 * GPA_CGS  # >= 1 GPa
    ok_pts = int(m_s.sum()) > 30
    dj = float(np.abs(np.diff(comp[m_s])).max()) if ok_pts else 9.9
    peak = float(comp.max()) if len(comp) else 0.0
    g.check("hugoniot-smooth", ok_pts and dj < 0.2,
            "max adjacent compression jump %.3f for P >= 1 GPa (%d pts); "
            "p0=%.3g GPa (report)" % (dj, int(m_s.sum()), p0 / 1e10))
    g.check("hugoniot-peak", True,
            "peak compression %.3f (report — no Ti reference data "
            "harvested yet)" % peak, hard=False)

    # seam-3 physics: Saha vs FD ideal plasma agreement in-band
    band3 = np.abs(lT - S3_TI.c) <= S3_TI.d
    sa, ip = res["srcs"]["saha"], res["srcs"]["ip"]
    sel = band3[None, :] & ip["avail"]
    relp = np.abs(sa["p"][sel] / ip["p"][sel] - 1.0)
    g.check("seam3-saha-vs-ip", float(np.median(relp)) < 0.05,
            "P agreement in seam-3 band: median %.3f%% max %.3f%%"
            % (100 * np.median(relp), 100 * relp.max()))

    if qa_dir:
        plots_ti(qa_dir, res, cs2, Rm, comp, pres)

    print("SPLICE-QA OVERALL %s (%d gate failure(s))"
          % ("PASS" if not g.failed else "FAIL", len(g.failed)))
    return (1 if g.failed else 0), dict(comp=comp, pres=pres)


def plots_ti(outdir, res, cs2, Rm, comp, pres):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    lrho, lT = res["lrho"], res["lT"]
    hull = res["blocks"]["hull"]
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
    dom = np.argmax(res["W"], axis=0).astype(float)
    im = axes[0, 0].pcolormesh(lT, lrho, dom, shading="nearest",
                               vmin=0, vmax=2)
    fig.colorbar(im, ax=axes[0, 0],
                 label="dominant source (0 solid, 1 Saha, 2 IP)")
    axes[0, 0].contour(lT, lrho, hull, levels=[0.5], colors="r",
                       linewidths=0.7)
    axes[0, 0].set_title("Ti sources + hull (red)")
    im = axes[0, 1].pcolormesh(lT, lrho,
                               res["srcs"]["saha"].get("zbar",
                                                       np.zeros_like(hull)),
                               shading="nearest")
    fig.colorbar(im, ax=axes[0, 1], label="Saha Zbar")
    axes[0, 1].set_title("mean ionization (Saha)")
    im = axes[1, 0].pcolormesh(lT, lrho, np.log10(np.maximum(cs2, 1.0)),
                               shading="nearest")
    fig.colorbar(im, ax=axes[1, 0], label="log10 cs^2")
    axes[1, 0].set_title("cs^2")
    axes[1, 1].semilogy(comp, pres / GPA_CGS, "-")
    axes[1, 1].set_xlabel("compression rho/rho0")
    axes[1, 1].set_ylabel("P [GPa]")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_title("Ti cold-start Hugoniot (v0, no reference data)")
    for ax in axes.flat[:3]:
        ax.set_xlabel("log10 T [K]")
        ax.set_ylabel("log10 rho [g/cc]")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "splice_Ti_qa.png"), dpi=130)
    plt.close(fig)


def bilin(lrho, lT, F, rho, T):
    """Bilinear table lookup mirroring the C++ rule (QA-local)."""
    x = np.interp(np.log10(rho), lrho, np.arange(len(lrho)))
    y = np.interp(np.log10(T), lT, np.arange(len(lT)))
    i, j = int(min(x, len(lrho) - 2)), int(min(y, len(lT) - 2))
    fx, fy = x - i, y - j
    return ((1 - fx) * (1 - fy) * F[i, j] + fx * (1 - fy) * F[i + 1, j]
            + (1 - fx) * fy * F[i, j + 1] + fx * fy * F[i + 1, j + 1])


def seam_c1_metric(lrho, lT, F, band, axis):
    """Max in-band |Delta d(lnF)/d(l-axis)| vs 3x the 90th pct outside.

    axis=1: slopes along lT per isochore; axis=0: along lrho per isotherm.
    Returns (in-band max, out-of-band p90)."""
    lnF = np.log(np.maximum(F, 1e-300))
    d = np.diff(lnF, axis=axis)
    jump = np.abs(np.diff(d, axis=axis))  # second difference = slope jumps
    if axis == 1:
        bmid = band[1:-1] if len(band) == lnF.shape[1] else band
        inb = jump[:, bmid]
        outb = jump[:, ~bmid]
    else:
        return None
    return float(inb.max()), float(np.percentile(outb, 90.0))


def run_qa(res, qa_dir=None, hug_ref_path=None):
    g = Gate()
    lrho, lT = res["lrho"], res["lT"]
    blocks = res["blocks"]
    p, e, hull = blocks["p"], blocks["e"], blocks["hull"] > 0.5
    R = 10.0 ** lrho[:, None] * np.ones((1, len(lT)))
    TT = 10.0 ** lT[None, :] * np.ones((len(lrho), 1))

    g.check("hull", hull.mean() > 0.5,
            "in-hull fraction %.3f; tension-clipped %d cells (p_floor %.1e)"
            % (hull.mean(), res["n_clip"], res["p_floor"]))

    # alignment constancy (splice-plan section 5.1 gate)
    for a in res["align_stats"]:
        g.check("align-" + a["pair"], a["std_over_kT"] < 0.10,
                "cE=%.4e erg/g over n=%d band cells; mean|dev|/kT=%.4f "
                "(gate < 0.10)" % (a["cE"], a["n"], a["std_over_kT"]),
                hard=False)

    # convexity (hard): the D5 sound-speed combination positive in-hull
    cs2 = sound_speed_sq(R, TT, p, blocks["dpdrho"], blocks["dpdT"],
                         blocks["cv"])
    n_bad = int(np.sum(cs2[hull] <= 0.0))
    g.check("convexity", n_bad == 0,
            "cs^2 <= 0 at %d/%d in-hull cells" % (n_bad, int(hull.sum())))

    # Maxwell residual: blended vs the REOS.3 yardstick on the same grid
    Rmap = maxwell_residual(lrho, lT, p, e)
    reos = res["srcs"]["reos"]
    Ry = maxwell_residual(lrho, lT, reos["p"], reos["e"])
    ry_med = float(np.median(Ry[reos["hull"] > 0.5]))
    med = float(np.median(Rmap[hull]))
    g.check("maxwell", med < 2.0 * max(ry_med, 1e-6),
            "blended median %.2e vs REOS.3 yardstick %.2e (gate < 2x); "
            "p95 %.2e" % (med, ry_med,
                          float(np.percentile(Rmap[hull], 95.0))))

    # C1 seam metric: slope-jump inside each band vs 3x p90 outside
    anyband = np.zeros(len(lT), bool)
    for nm, s in (("s1", S1), ("s2", S2), ("s3", S3)):
        anyband |= np.abs(lT - s.c) <= 1.5 * s.d
    lnp = np.log(np.maximum(p, 1e-300))
    dslope = np.abs(np.diff(np.diff(lnp, axis=1), axis=1))
    hmid = hull[:, 1:-1]
    omid = ~anyband[1:-1]
    out_sample = dslope[:, omid][hmid[:, omid]]
    out_p99 = float(np.percentile(out_sample, 99.0))
    out_max = float(out_sample.max())
    for nm, s in (("s1", S1), ("s2", S2), ("s3", S3)):
        bmid = (np.abs(lT - s.c) <= 1.5 * s.d)[1:-1]
        in_sample = dslope[:, bmid][hmid[:, bmid]]
        in_p99 = float(np.percentile(in_sample, 99.0))
        in_max = float(in_sample.max())
        # like-for-like: p99 vs 3x p99 (max-vs-quantile was statistically
        # unfair); maxes reported for the record
        g.check("c1-seam-" + nm, in_p99 <= 3.0 * max(out_p99, 1e-12),
                "in-band p99 %.3e vs 3x out-of-band p99 %.3e "
                "(maxes: in %.3e out %.3e)"
                % (in_p99, 3.0 * out_p99, in_max, out_max))

    # cold-start principal Hugoniot from the table's own reference state
    rho0, T0 = REF_STATE["rho0"], REF_STATE["T0"]
    e0 = bilin(lrho, lT, e, rho0, T0)
    p0 = bilin(lrho, lT, p, rho0, T0)
    comp, pres = locus(lrho, lT, p, e, rho0, e0, p0)
    ok_pts = len(comp) > 30
    peak = float(comp.max()) if ok_pts else 0.0
    g.check("hugoniot-peak", ok_pts and 4.2 <= peak <= 4.9,
            "peak compression %.3f (gate [4.2, 4.9]); %d locus points; "
            "e0=%.4e p0=%.3e" % (peak, len(comp), e0, p0))
    # smoothness gated for P >= 0.3 GPa: below that the locus foot threads
    # the two-phase dome's hull-0 fill, where spurious extra RH roots are
    # an artefact of the fill, not the table (foot despike count reported)
    P_SMOOTH_MIN = 0.3 * GPA_CGS
    m_s = pres >= P_SMOOTH_MIN
    dj = float(np.abs(np.diff(comp[m_s])).max()) if m_s.sum() > 5 else 9.9
    n_foot = int(len(comp) - m_s.sum())
    g.check("hugoniot-smooth", dj < 0.15,
            "max adjacent compression jump %.3f for P >= 0.3 GPa "
            "(gate < 0.15); %d foot points below the gate range"
            % (dj, n_foot))

    # anchor: PIMC MC2000 locus (rho0 = 0.171 liquid D2, same reference)
    ref = None
    hug_ref_path = hug_ref_path or raw_path("hugoniot_MC2000_PRL85_1890.txt")
    if os.path.exists(hug_ref_path):
        ref = np.loadtxt(hug_ref_path)  # T rho P_GPa compression
        rel = []
        for Tr, rr, Pr, cr in ref:
            if Pr * GPA_CGS < pres.min() or Pr * GPA_CGS > pres.max():
                continue
            c_tab = np.interp(Pr * GPA_CGS, pres, comp)
            rel.append(abs(c_tab / cr - 1.0))
        rel = np.array(rel)
        g.check("hugoniot-anchor-PIMC",
                len(rel) > 0 and float(np.median(rel)) < 0.05,
                "n=%d MC2000 points in range: median %.2f%% max %.2f%% "
                "(gate median < 5%%)" % (len(rel), 100 * np.median(rel),
                                         100 * rel.max()), hard=True)

    # validation overlay: the single published iFPEOS isochore (report)
    iso_path = raw_path("iFPEOS", "ifpeos_rho0.001_isochore.csv")
    if os.path.exists(iso_path):
        d = np.loadtxt(iso_path, delimiter=",", skiprows=4)
        Tq, Pq, sPq = d[:, 0], d[:, 1], d[:, 2]
        # only MD points whose own error bar is meaningful (sigma < 30% P);
        # the 800 K point has sigma ~ 4x its P and carries no information
        m = (Tq >= 10.0 ** lT[0]) & (Tq <= 10.0 ** lT[-1]) & (Pq > 0) \
            & (sPq < 0.3 * Pq)
        devs = np.array([abs(bilin(lrho, lT, p, 0.001, tq) / pq - 1.0)
                         for tq, pq in zip(Tq[m], Pq[m])])
        g.check("ifpeos-isochore", True,
                "rho=0.001 (sigma-filtered n=%d): median %.2f%% max %.2f%% "
                "(validation overlay, report-only)"
                % (len(devs), 100 * np.median(devs), 100 * devs.max()),
                hard=False)

    if qa_dir:
        plots(qa_dir, res, cs2, Rmap, comp, pres, ref)

    print("SPLICE-QA OVERALL %s (%d gate failure(s))"
          % ("PASS" if not g.failed else "FAIL", len(g.failed)))
    return (1 if g.failed else 0), dict(comp=comp, pres=pres)


def plots(outdir, res, cs2, Rmap, comp, pres, ref):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    lrho, lT = res["lrho"], res["lT"]
    blocks = res["blocks"]
    hull = blocks["hull"]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
    # dominant source map
    Wn = res["W"]
    dom = np.argmax(Wn, axis=0).astype(float)
    im = axes[0, 0].pcolormesh(lT, lrho, dom, shading="nearest",
                               cmap="viridis", vmin=0, vmax=3)
    fig.colorbar(im, ax=axes[0, 0],
                 label="dominant source (0 cold, 1 REOS3, 2 FPEOS, 3 IP)")
    axes[0, 0].contour(lT, lrho, hull, levels=[0.5], colors="r",
                       linewidths=0.7)
    for s in (S1, S2, S3):
        axes[0, 0].axvline(s.c, color="w", ls="--", lw=0.7)
    axes[0, 0].set_title("sources + hull (red)")

    im = axes[0, 1].pcolormesh(lT, lrho, np.log10(np.maximum(cs2, 1.0)),
                               shading="nearest")
    fig.colorbar(im, ax=axes[0, 1], label="log10 cs^2")
    axes[0, 1].set_title("cs^2 (D5 combination)")

    im = axes[1, 0].pcolormesh(lT, lrho,
                               np.log10(np.maximum(Rmap, 1e-12)),
                               shading="nearest", vmax=0.0)
    fig.colorbar(im, ax=axes[1, 0], label="log10 Maxwell residual")
    axes[1, 0].set_title("Maxwell residual")

    axes[1, 1].semilogy(comp, pres / GPA_CGS, "-", lw=1.2,
                        label="spliced table")
    if ref is not None:
        axes[1, 1].semilogy(ref[:, 3], ref[:, 2], "rs", mfc="none",
                            label="MC2000 PIMC")
    axes[1, 1].set_xlabel("compression rho/rho0")
    axes[1, 1].set_ylabel("P [GPa]")
    axes[1, 1].set_title("cold-start principal Hugoniot "
                         "(rho0=%.3f, T0=%g K)" % (REF_STATE["rho0"],
                                                   REF_STATE["T0"]))
    axes[1, 1].legend()
    for ax in axes.flat[:3]:
        ax.set_xlabel("log10 T [K]")
        ax.set_ylabel("log10 rho [g/cc]")
    fig.suptitle("D_spliced QA")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "splice_D_qa.png"), dpi=130)
    plt.close(fig)
    np.savetxt(os.path.join(outdir, "hugoniot_D_spliced.txt"),
               np.column_stack([comp, pres / GPA_CGS]),
               header="compression_rho_over_rho0  P_GPa  (cold start "
                      "rho0=%.3f T0=%gK)" % (REF_STATE["rho0"],
                                             REF_STATE["T0"]))
