"""QA plots and printed diagnostics for conditioned tables."""

import os

import numpy as np
from scipy.interpolate import PchipInterpolator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .constants import GPA_CGS
from .hugoniot import locus


def qa_plots(outdir, lrho, lT, p, e, cv, hull, pts=None, hug=None, loo=None,
             tag=""):
    os.makedirs(outdir, exist_ok=True)
    rho = 10.0 ** lrho
    T = 10.0 ** lT

    # isotherms P(rho), e(rho) with raw points overplotted
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    jsel = np.unique(np.linspace(0, len(lT) - 1, 10).astype(int))
    for j in jsel:
        axes[0].loglog(rho, p[:, j], label="T=%.3g K" % T[j])
        axes[1].loglog(rho, e[:, j])
    if pts is not None:
        for j in jsel:
            m = np.abs(np.log10(pts["T"]) - lT[j]) < 0.02
            axes[0].plot(pts["rho"][m], pts["p"][m], "k.", ms=4)
            axes[1].plot(pts["rho"][m], pts["e"][m] - pts["e"].min() + 1e-3
                         * (pts["e"].max() - pts["e"].min()), "k.", ms=4)
    axes[0].set_xlabel("rho [g/cc]")
    axes[0].set_ylabel("P [erg/cc]")
    axes[1].set_xlabel("rho [g/cc]")
    axes[1].set_ylabel("e (shifted) [erg/g]")
    axes[0].legend(fontsize=6)
    fig.suptitle("isotherms " + tag)
    fig.savefig(os.path.join(outdir, "isotherms%s.png" % tag), dpi=130)
    plt.close(fig)

    # cv heatmap + hull
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    im = axes[0].pcolormesh(lT, lrho, np.log10(cv), shading="nearest")
    fig.colorbar(im, ax=axes[0], label="log10 cv [erg/g/K]")
    axes[1].pcolormesh(lT, lrho, hull, shading="nearest", cmap="gray")
    for ax, t in zip(axes, ("cv", "hull mask (white=inside)")):
        ax.set_xlabel("log10 T [K]")
        ax.set_ylabel("log10 rho [g/cc]")
        ax.set_title(t)
    fig.savefig(os.path.join(outdir, "cv_hull%s.png" % tag), dpi=130)
    plt.close(fig)

    # round-trip inversion residuals (python prototype of Stage-2 C++)
    ee = np.empty_like(e)
    for i in range(len(lrho)):
        fe = PchipInterpolator(lT, e[i])
        for j in range(len(lT)):
            # bisection solve e(T)=e[i,j] (monotone segments assumed piecewise)
            lo, hi = lT[0], lT[-1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if fe(mid) < e[i, j]:
                    lo = mid
                else:
                    hi = mid
            ee[i, j] = fe(0.5 * (lo + hi))
    res = np.abs(ee - e) / np.maximum(np.abs(e), 1e-30)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.log10(np.maximum(res[hull > 0.5], 1e-18)), bins=60)
    ax.set_xlabel("log10 relative e-residual (rt->re->rt prototype)")
    ax.set_title("max=%.2e  (hull cells)" % res[hull > 0.5].max())
    fig.savefig(os.path.join(outdir, "roundtrip%s.png" % tag), dpi=130)
    plt.close(fig)

    # Hugoniot from the conditioned table
    if hug is not None:
        rho0, e0, p0 = hug["rho0"], hug["e0"], hug["p0"]
        hr, hp = locus(lrho, lT, p, e, rho0, e0, p0)
        hr, hp = list(hr), list(hp)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.semilogy(hr, np.array(hp) / GPA_CGS, "o", ms=3, label="this table")
        if hug.get("ref") is not None:
            ref = hug["ref"]  # columns: T rho P_GPa compression
            ax.semilogy(ref[:, 3], ref[:, 2], "rs", mfc="none",
                        label=hug.get("ref_label", "published"))
        ax.set_xlabel("compression rho/rho0")
        ax.set_ylabel("P [GPa]")
        ax.set_title("principal Hugoniot (rho0=%.4g g/cc)" % rho0)
        ax.legend()
        fig.savefig(os.path.join(outdir, "hugoniot%s.png" % tag), dpi=130)
        plt.close(fig)
        np.savetxt(os.path.join(outdir, "hugoniot%s.txt" % tag),
                   np.column_stack([hr, np.array(hp) / GPA_CGS]),
                   header="compression_rho_over_rho0  P_GPa")

    if loo is not None and len(loo):
        print("leave-one-out regrid error (isochore pass): "
              "P median %.3g max %.3g | e median %.3g max %.3g"
              % (np.median(loo[:, 0]), loo[:, 0].max(),
                 np.median(loo[:, 1]), loo[:, 1].max()))
