"""SS2 QA for the deuterium cold/low-T composite (plan §3, SS2 gates).

Emits one 'COLDMODEL-QA[<gate>] PASS/FAIL <numbers>' line per gate (the
same greppable convention as the C++ self-test) plus plots. Hard gates set
the exit code; 'report' lines are informational.
"""

import io
import os
import tarfile

import numpy as np

from .condition import condition
from .constants import GPA_CGS, KB, M_D, NA  # noqa: F401 (NA used in Ti QA)
from .formats.fpeos import fpeos_to_deuterium, read_fpeos
from .materials.deuterium import RHO0, T0, DeuteriumColdModel
from .models.ideal_plasma import IdealPlasma
from .thermo import maxwell_residual, sound_speed_sq
from . import sources as _sources

# QA evaluation rectangle (log10): the cold model owns lT < ~3.4 (seam-1
# top + margin); rho spans the table range up to past the solid ramp
LRHO = (-4.0, 0.8)
LT = (1.25, 3.6)


class Gate:
    def __init__(self):
        self.failed = []

    def check(self, name, ok, detail, hard=True):
        verdict = "PASS" if ok else ("FAIL" if hard else "FAIL(report)")
        print("COLDMODEL-QA[%s] %s %s" % (name, verdict, detail))
        if hard and not ok:
            self.failed.append(name)


def fpeos_points_dbasis():
    """FPEOS H table read straight out of the committed tarball, D basis."""
    tgz = os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                       "data", "raw", "fpeos_10-26-25.tar.gz")
    with tarfile.open(tgz, "r:gz") as tf:
        raw = tf.extractfile("FPEOS/H_EOS_09-18-20.txt").read().decode()
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(raw)
        path = f.name
    try:
        return fpeos_to_deuterium(read_fpeos(path))
    finally:
        os.unlink(path)


def run_qa(args):
    if args.material == "Ti":
        return run_qa_ti(args)
    g = Gate()
    model = DeuteriumColdModel(use_coolprop=not args.no_coolprop)

    # --- gate: reference-state closure (CoolProp = Richardson 2014) ------
    if model.cp is not None:
        import CoolProp.CoolProp as CP
        rho_ref = CP.PropsSI("D", "T", T0, "P", 1.0e5, "Deuterium") * 1e-3
        cs_ref = CP.PropsSI("A", "T", T0, "P", 1.0e5, "Deuterium") * 1e2
        g.check("ref-density", abs(rho_ref / RHO0 - 1.0) < 5e-3,
                "rho(20K,1bar)=%.5f g/cc vs %.3f (tol 0.5%%)" % (rho_ref, RHO0))
        g.check("ref-sound-speed", 0.9 * 1.1e5 < cs_ref < 1.1 * 1.1e5,
                "cs(20K,1bar)=%.4g cm/s vs ~1.1 km/s (tol 10%%)" % cs_ref)

    # --- gate: cold-curve fit quality ------------------------------------
    g.check("vinet-fit", model.vinet_rms < 0.15,
            "rms rel P err %.3f over LLNL 0-100 GPa (v0=%.4g cc/g, "
            "B0=%.4g GPa, B0'=%.3f, rho_switch=%.3f g/cc)"
            % (model.vinet_rms, model.vinet[0], model.vinet[1] / GPA_CGS,
               model.vinet[2], model.rho_switch))

    # --- gate: alignment constancy ---------------------------------------
    for key, st in model.align_stats.items():
        g.check("align-" + key,
                st["std_e_over_kT"] < 0.5 and st["std_s_over_R"] < 0.5,
                "cE=%.3e erg/g cS=%.3e erg/g/K; std_e/kT=%.3f std_s/R=%.3f"
                % (st["cE"], st["cS"], st["std_e_over_kT"],
                   st["std_s_over_R"]), hard=False)

    # --- grid evaluation --------------------------------------------------
    lrho = np.linspace(*LRHO, args.n_rho)
    lT = np.linspace(*LT, args.n_T)
    R, T = np.meshgrid(10.0 ** lrho, 10.0 ** lT, indexing="ij")
    out = model.eval(R, T)
    hull = out["hull"] > 0.5
    finite = all(np.all(np.isfinite(out[k])) for k in ("p", "e", "s", "a"))
    g.check("finite", finite, "all p/e/s/a finite on %dx%d grid (hull %.1f%%)"
            % (args.n_rho, args.n_T, 100.0 * hull.mean()))

    # --- gates: convexity + cv > 0 via the production conditioning -------
    e_off = out["e"] - out["e"][hull].min() + 0.05 * np.ptp(out["e"][hull])
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, out["p"], e_off)
    cs2 = sound_speed_sq(R, T, out["p"], dpdrho, dpdT, cv)
    n_cs = int(np.sum(cs2[hull] <= 0.0))
    g.check("convexity", n_cs == 0,
            "cs^2 <= 0 at %d/%d in-hull cells" % (n_cs, int(hull.sum())))
    cv_raw_neg = stats["cv_floored"]
    g.check("cv-positive", True,
            "cv floored at %d cells (incl. off-hull; floor %.3g)"
            % (cv_raw_neg, stats["cv_floor"]), hard=False)
    n_cv_hull = int(np.sum(cv[hull] <= 0.0))
    g.check("cv-positive-hull", n_cv_hull == 0,
            "cv <= 0 at %d in-hull cells (post-floor)" % n_cv_hull)

    # --- gate: Maxwell residual (analytic composite) ----------------------
    Rmap = maxwell_residual(lrho, lT, out["p"], out["e"])
    med = float(np.median(Rmap[hull]))
    p95 = float(np.percentile(Rmap[hull], 95.0))
    g.check("maxwell", med < 1e-3,
            "median %.2e p95 %.2e in-hull (gate: median < 1e-3; p95 "
            "reported — includes grid-PCHIP truncation in the metric)"
            % (med, p95))

    # --- gate: ideal-plasma anchor vs FPEOS at T >= 3e7 K ------------------
    pts = fpeos_points_dbasis()
    m = pts["T"] >= 3.0e7
    ip = IdealPlasma(M_D, Z=1.0)
    p_mod = ip.p(pts["rho"][m], pts["T"][m])
    rel = np.abs(p_mod / pts["p"][m] - 1.0)
    g.check("anchor-vs-fpeos",
            float(np.median(rel)) < 0.02 and float(rel.max()) < 0.05,
            "n=%d T>=3e7K: median %.3f%% max %.3f%% (gate median<2%%, "
            "max<5%%)" % (int(m.sum()), 100 * np.median(rel), 100 * rel.max()))

    # --- plots -------------------------------------------------------------
    if args.qa:
        plots(args.qa, lrho, lT, out, cv, cs2, Rmap, model)

    print("COLDMODEL-QA OVERALL %s (%d gate failure(s))"
          % ("PASS" if not g.failed else "FAIL", len(g.failed)))
    return 1 if g.failed else 0


def run_qa_ti(args):
    """SS5b Ti solid-model QA (SS2-prime): the same gate classes as the
    deuterium cold model, on the harvested LLNL Ti 0 K isotherm.
    Literature values quoted in reports (theta_D ~ 420 K, B0 ~ 110 GPa,
    gamma_e ~ 3.3-3.6 mJ/mol/K^2) are cross-checks, never inputs."""
    import numpy as np

    from .materials.titanium import M_TI, RHO0_AMBIENT, TitaniumColdModel

    g = Gate()
    model = TitaniumColdModel()
    v0, B0, B0p = model.vinet

    g.check("vinet-fit", model.vinet_rms < 0.05,
            "rms rel P err %.4f over LLNL Ti 0-100 GPa (v0=%.4f cc/g, "
            "B0=%.1f GPa, B0'=%.3f)" % (model.vinet_rms, v0,
                                        B0 / GPA_CGS, B0p))
    rho0_fit = 1.0 / v0
    g.check("rho0-closure", abs(rho0_fit / 4.580 - 1.0) < 0.01,
            "rho0(fit)=%.4f g/cc vs 4.580 (LLNL 0K datum; ambient 300K "
            "%.3f)" % (rho0_fit, RHO0_AMBIENT))
    th0 = model.theta0()
    g.check("theta0-slater", True,
            "Slater theta0 = %.0f K (literature Debye ~420 K; "
            "report-only — Slater ignores shear)" % th0, hard=False)
    ge = float(model.electrons.gamma_spec(rho0_fit)) * M_TI * NA * 1e-7 * 1e3
    g.check("gamma-e", True,
            "free-electron gamma_e(z_c=4) = %.2f mJ/mol/K^2 (literature "
            "~3.3-3.6 incl. d-band enhancement; report-only)" % ge,
            hard=False)

    # grid gates: solid Ti rectangle, ambient to ~5.5x compression
    lrho = np.linspace(np.log10(3.0), np.log10(25.0), args.n_rho)
    lT = np.linspace(1.25, 3.6, args.n_T)
    R, T = np.meshgrid(10.0 ** lrho, 10.0 ** lT, indexing="ij")
    out = model.eval(R, T)
    finite = all(np.all(np.isfinite(out[k])) for k in ("p", "e", "s", "a"))
    g.check("finite", finite, "all p/e/s/a finite on %dx%d grid"
            % (args.n_rho, args.n_T))

    e_off = out["e"] - out["e"].min() + 0.05 * np.ptp(out["e"])
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, out["p"], e_off)
    cs2 = sound_speed_sq(R, T, out["p"], dpdrho, dpdT, cv)
    # tension region (p < 0 below rho0 at low T) is excluded: the splice
    # clips it; the physical gate is the compressed branch
    phys = out["p"] > 0.0
    n_cs = int(np.sum(cs2[phys] <= 0.0))
    g.check("convexity", n_cs == 0,
            "cs^2 <= 0 at %d/%d p>0 cells" % (n_cs, int(phys.sum())))
    cs0 = np.sqrt(float(np.interp(np.log10(4.51),
                                  lrho, cs2[:, 0]))) / 1e5
    g.check("cs-ambient", True,
            "cs(4.51 g/cc, ~18 K) = %.2f km/s (literature Ti bulk sound "
            "speed ~4.9-5.2; report-only)" % cs0, hard=False)

    Rmap = maxwell_residual(lrho, lT, out["p"], out["e"])
    med = float(np.median(Rmap[phys]))
    g.check("maxwell", med < 1e-3,
            "median %.2e p95 %.2e over p>0 cells (analytic model)"
            % (med, float(np.percentile(Rmap[phys], 95.0))))

    print("COLDMODEL-QA OVERALL %s (%d gate failure(s))"
          % ("PASS" if not g.failed else "FAIL", len(g.failed)))
    return 1 if g.failed else 0


def plots(outdir, lrho, lT, out, cv, cs2, Rmap, model):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    rho = 10.0 ** lrho
    T = 10.0 ** lT
    hull = out["hull"] > 0.5

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    jsel = np.unique(np.linspace(0, len(lT) - 1, 9).astype(int))
    for j in jsel:
        axes[0, 0].loglog(rho, np.maximum(out["p"][:, j], 1e2),
                          label="T=%.3g K" % T[j])
    # LLNL cold-curve overlay
    P_GPa, V_mol = np.loadtxt(
        os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                     "data", "raw", "llnl_coldcurve", "h_coldcurve_0K.csv"),
        delimiter=",", skiprows=4, unpack=True)
    rho_cc = 1.0 / (V_mol / (M_D * NA))
    axes[0, 0].loglog(rho_cc, np.maximum(P_GPa, 1e-8) * GPA_CGS, "k.", ms=4,
                      label="LLNL 0 K")
    axes[0, 0].set_xlabel("rho [g/cc]")
    axes[0, 0].set_ylabel("P [erg/cc]")
    axes[0, 0].legend(fontsize=6)
    axes[0, 0].set_title("isotherms + LLNL cold curve")

    im = axes[0, 1].pcolormesh(lT, lrho, out["src"] + 0.0, shading="nearest",
                               cmap="viridis")
    fig.colorbar(im, ax=axes[0, 1], label="src (0 CP, 1 gas, 2 solid)")
    axes[0, 1].contour(lT, lrho, out["hull"], levels=[0.5], colors="r")
    axes[0, 1].set_title("dominant source + hull (red)")

    im = axes[1, 0].pcolormesh(lT, lrho, np.log10(np.maximum(cs2, 1e0)),
                               shading="nearest")
    fig.colorbar(im, ax=axes[1, 0], label="log10 cs^2 [cgs]")
    axes[1, 0].set_title("sound speed squared")

    im = axes[1, 1].pcolormesh(
        lT, lrho, np.log10(np.maximum(Rmap, 1e-12)), shading="nearest",
        vmax=0.0)
    fig.colorbar(im, ax=axes[1, 1], label="log10 Maxwell residual")
    axes[1, 1].set_title("Maxwell/Grueneisen residual (in-hull med %.1e)"
                         % np.median(Rmap[hull]))
    for ax in axes.flat[1:]:
        ax.set_xlabel("log10 T [K]")
        ax.set_ylabel("log10 rho [g/cc]")
    fig.suptitle("D cold composite QA (Vinet rms %.3f, rho_switch %.3f g/cc)"
                 % (model.vinet_rms, model.rho_switch))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "coldmodel_D_qa.png"), dpi=130)
    plt.close(fig)
