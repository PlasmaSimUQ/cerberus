"""Cold extension of a SESAME backbone with an analytic solid model
(doc/ti_splice_plan_v2.md T3/T4, decisions B2/B3/B6).

Replaces the constant-in-T nearest-fill below the SESAME hull with a
physically-ordered construction, per target isochore (column i):

  zone 1  T <= min(F*T_m, T_hull_lo):  solid model values, hull 0
          (model-supplied, not source-supplied);
  zone 2  compressed columns (T_hull_lo < F*T_m): tanh blend in log T
          between solid and SESAME over >= 2*N_HB cells (R3), hull kept;
          the H5 energy offset is measured here;
  zone 3  expanded columns (T_hull_lo > T_m, e.g. the 4.40 g/cc casing):
          linear-in-logT bridge between the solid anchor at F*T_m and the
          SESAME anchor at T_hull_lo, hull 0 — the declared model-free
          bridge across SESAME's own liquid gap.

The melt cap does double duty: T_m -> 0 at low density, so the solid
model is automatically excluded from the vapour region (no repeat of the
v0 ideal-vapour fill, R1). Columns where F*T_m falls below the grid's T
floor are left untouched (trackP semantics unchanged there).

A second crossover pass (maxwell.crossover_isotherm) then repairs the
solid model's tension foot (p < 0 at rho < rho0) with the standard
monotone ramp + lever-rule e; those cells join the band (-> c_cav floor).

Energy alignment (H5): ONE constant offset solid->SESAME, measured over
the zone-2 overlap; gated on constancy (std over the local kT scale).
Pressures are never shifted; the p mismatch in the overlap is reported.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

from .constants import GPA_CGS, KB

N_HB = 10          # blend half-band in cells (R3: >= 10 per half-band)
W_EPS = 1e-4       # weights below this are zeroed (as splice.py)


def melt_curve_on(lrho, ml):
    """T_m(rho) from a 411 dict (read_sesame_1d), PCHIP in log10 rho,
    clamped to the 411's own positive-density range."""
    rho = np.asarray(ml["rho"], float)
    Tm = np.asarray(ml["Tm"], float)
    m = rho > 0.0
    lr = np.log10(rho[m])
    f = PchipInterpolator(lr, Tm[m])
    x = np.clip(np.asarray(lrho, float), lr[0], lr[-1])
    return np.asarray(f(x), float)


def check_cold_model(model, rho0, cc, p_tol_gpa=0.5):
    """S4 acceptance gates on the fitted solid model (plan v2).

    - |p(rho0, 300 K)| <= p_tol_gpa (the 0 K density sits just above rho0,
      so the thermal pressure must close the 306's tension there);
    - fitted 0 K density within 10% of the 201 rho0;
    - B0 within 30% of the raw 306 slope stiffness at rho0.
    Returns a report dict; raises RuntimeError on a hard failure.
    """
    v0, B0, B0p = model.vinet
    p_amb = float(model.p(rho0, 300.0)) / GPA_CGS
    rho = np.asarray(cc["rho"], float)
    P = np.asarray(cc["p"], float)
    i = int(np.argmin(np.abs(rho - rho0)))
    i = min(max(i, 1), len(rho) - 2)
    B306 = rho0 * (P[i + 1] - P[i - 1]) / (rho[i + 1] - rho[i - 1])  # GPa
    rep = dict(rho0K=1.0 / v0, B0_gpa=B0 / GPA_CGS, B0p=B0p,
               p_amb_gpa=p_amb, B306_gpa=B306, theta0=model.theta0())
    fails = []
    if abs(p_amb) > p_tol_gpa:
        fails.append("p(rho0,300K) = %+.3f GPa exceeds +-%.2f" % (p_amb, p_tol_gpa))
    if abs(rep["rho0K"] / rho0 - 1.0) > 0.10:
        fails.append("fitted 0K density %.3f vs 201 rho0 %.3f" % (rep["rho0K"], rho0))
    if not (0.7 <= rep["B0_gpa"] / max(B306, 1e-30) <= 1.3):
        fails.append("B0 %.1f GPa vs raw-306 slope %.1f GPa" % (rep["B0_gpa"], B306))
    if fails:
        raise RuntimeError("cold-model gate FAILED: " + "; ".join(fails))
    return rep


def cold_extend(lrho, lT, p, e, hull, band, model, Tm, f_melt=0.9,
                align_tol=5.0, m_ref=None, rho_max=None):
    """Apply the three-zone extension in place on copies; returns
    (p, e, hull, band, stats). Arrays are (Nr, Nt); hull/band as in
    cmd_sesame after the band remap. Tm is T_m on the lrho axis.

    rho_max: solid-model validity ceiling in g/cc — normally the Vinet fit
    window's upper edge. MEASURED on 2963 (plan v2 T3): inside the fit
    window the solid<->SESAME energy offset is constant to ~3.5 kT across
    columns (per-column std/kT <= 0.14, p mismatch <= 9%); beyond it the
    offset drifts (1.2e12 erg/g by 50 g/cc, p mismatch 64%) — blending an
    extrapolated model there would corrupt cv in the blend band. Columns
    above rho_max keep the backbone fill, which below ~75 K is close to
    physical anyway (a cold solid has dp/dT ~ 0 for T << theta)."""
    p = p.copy()
    e = e.copy()
    hull = hull.copy()
    band = band.copy()
    nr, nt = p.shape
    T = 10.0 ** np.asarray(lT, float)
    dlT = lT[1] - lT[0]

    # per-column indices
    j_hull = np.full(nr, -1, int)
    for i in range(nr):
        ins = np.flatnonzero(hull[i] > 0.5)
        if ins.size:
            j_hull[i] = ins[0]
    # last index with T <= f*Tm (-1 if none)
    j_cap = np.searchsorted(T, f_melt * np.asarray(Tm, float), side="right") - 1

    ext = (j_cap >= 0) & (j_hull >= 0)
    if rho_max is not None:
        ext &= (10.0 ** np.asarray(lrho, float)) <= rho_max
    if not ext.any():
        return p, e, hull, band, dict(cols=0)

    # solid-model evaluation on the needed subrectangle only (the Debye
    # integral is per-cell quadrature — keep it to the cold corner).
    # Overlap columns need the blend to decay below W_EPS (tanh arg ~4.6
    # half-widths past the centre -> jh + 7*N_HB is ample); gap columns
    # need the cap cell. Both are bounded by min(j_cap, jh + 7*N_HB).
    jmax = int(min(nt - 1,
                   np.max(np.minimum(j_cap[ext], j_hull[ext] + 7 * N_HB))))
    rows = np.flatnonzero(ext)
    R = (10.0 ** lrho[rows])[:, None] * np.ones((1, jmax + 1))
    TT = T[None, :jmax + 1] * np.ones((len(rows), 1))
    ps = np.asarray(model.p(R, TT), float)
    es = np.asarray(model.e(R, TT), float)

    # --- H5 energy alignment over the zone-2 overlap --------------------
    de, kts = [], []
    m_ref = m_ref or model.m_atom
    for k, i in enumerate(rows):
        jh = j_hull[i]
        if jh > j_cap[i] or jh < 0:      # gap column — no overlap here
            continue
        j1 = min(jh + 2 * N_HB, j_cap[i], jmax)
        for j in range(jh, j1 + 1):
            de.append(e[i, j] - es[k, j])
            kts.append(1.5 * KB * T[j] / m_ref)
    stats = dict(cols=int(ext.sum()))
    if de:
        de = np.asarray(de)
        cE = float(de.mean())
        std_over_kT = float(np.mean(np.abs(de - cE) / np.asarray(kts)))
        stats.update(align_n=len(de), align_cE=cE, align_std_kT=std_over_kT)
        if std_over_kT > align_tol:
            raise RuntimeError(
                "H5 offset-constancy gate FAILED: std/kT = %.3f > %.2f "
                "(n=%d, cE=%.4e erg/g) — solid model and SESAME disagree "
                "beyond a constant shift in the overlap" %
                (std_over_kT, align_tol, len(de), cE))
        es = es + cE
    else:
        stats.update(align_n=0, align_cE=0.0, align_std_kT=0.0)

    # p is never shifted; report the overlap mismatch instead (QA)
    pmis = []
    for k, i in enumerate(rows):
        jh = j_hull[i]
        if jh > j_cap[i] or jh < 0:
            continue
        j1 = min(jh + 2 * N_HB, j_cap[i], jmax)
        sl = slice(jh, j1 + 1)
        pmis.append(np.max(np.abs(ps[k, sl] - p[i, sl])
                           / np.maximum(np.abs(p[i, sl]), 1e-300)))
    stats["p_mismatch_max"] = float(max(pmis)) if pmis else 0.0

    # --- apply zones per column -----------------------------------------
    n_z1 = n_bl = n_br = n_trunc = 0
    for k, i in enumerate(rows):
        jh, jc = j_hull[i], int(j_cap[i])
        if jh <= jc:
            # overlap column: pure solid below the hull edge, tanh blend
            # centred N_HB cells above it, applied until the weight decays
            # below W_EPS naturally (no truncation seam; R3 half-band
            # = N_HB cells). The melt cap may cut it short (counted).
            c = lT[jh] + N_HB * dlT
            hw = N_HB * dlT
            jtop = min(jc, jmax)
            for j in range(0, jtop + 1):
                w = 0.5 * (1.0 - np.tanh((lT[j] - c) / hw))
                if j < jh:
                    w = 1.0
                elif w < W_EPS:
                    break                # tanh fully decayed — clean exit
                p[i, j] = w * ps[k, j] + (1.0 - w) * p[i, j]
                e[i, j] = w * es[k, j] + (1.0 - w) * e[i, j]
                if j < jh:
                    hull[i, j] = 0.0     # model-supplied
                    n_z1 += 1
                else:
                    n_bl += 1            # blend keeps source hull
            else:
                # loop ran to the melt cap with weight still > W_EPS
                n_trunc += 1
        else:
            # gap column: solid to the cap, declared bridge to the hull edge
            for j in range(0, jc + 1):
                p[i, j] = ps[k, j]
                e[i, j] = es[k, j]
                hull[i, j] = 0.0
                n_z1 += 1
            t = (lT[jc + 1:jh] - lT[jc]) / (lT[jh] - lT[jc])
            p[i, jc + 1:jh] = ps[k, jc] + t * (p[i, jh] - ps[k, jc])
            e[i, jc + 1:jh] = es[k, jc] + t * (e[i, jh] - es[k, jc])
            hull[i, jc + 1:jh] = 0.0
            n_br += int(jh - jc - 1)
    stats.update(zone1=n_z1, blend=n_bl, bridge=n_br, blend_trunc=n_trunc)
    return p, e, hull, band, stats


def cold_extend_stage(src, mat_id, raw, lrho, lT, p, e, hull, band,
                      f_melt=0.9, align_tol=5.0):
    """Full T3/T4 stage for cmd_sesame: read 306/411, fit + gate the solid
    model, extend, then the second crossover pass (plan v2 §0.3).

    Hull policy on the second pass: cells the crossover moves by more than
    1e-6 relative (or that were nonpositive) are banded/hull-0; smaller
    moves are regrid-ripple repair and keep their hull — the same policy
    the post-regrid monotonisation applies.
    """
    from .constants import AMU_G
    from .formats.sesame_ascii2 import read_sesame_1d
    from .materials.titanium import TitaniumColdModel, fit_from_306
    from .maxwell import condition_surface

    if not raw["abar"] or not raw["rho0"]:
        raise RuntimeError("--cold-extend needs the material's 201 table "
                           "(abar, rho0)")
    cc = read_sesame_1d(src, mat_id, 306)
    ml = read_sesame_1d(src, mat_id, 411)
    fit_rho = (4.0, 12.0)             # Vinet window = extension validity
    vinet, rms, n_fit = fit_from_306(cc, *fit_rho)
    model = TitaniumColdModel(vinet=vinet, m_atom=float(raw["abar"]) * AMU_G,
                              vinet_rms=rms)
    rep = check_cold_model(model, float(raw["rho0"]), cc)
    print("cold model (306 fit, n=%d, rms=%.2e): rho0K=%.4f g/cc  "
          "B0=%.1f GPa (raw306 %.1f)  B0p=%.3f  p(rho0,300K)=%+.3f GPa  "
          "theta0=%.0f K"
          % (n_fit, rms, rep["rho0K"], rep["B0_gpa"], rep["B306_gpa"],
             rep["B0p"], rep["p_amb_gpa"], rep["theta0"]))

    Tm = melt_curve_on(lrho, ml)
    p, e, hull, band, st = cold_extend(lrho, lT, p, e, hull, band, model,
                                       Tm, f_melt=f_melt,
                                       align_tol=align_tol,
                                       rho_max=fit_rho[1])
    print("cold-extend (melt cap 411, f=%.2f, rho<=%g): %d cols; zone1=%d "
          "blend=%d bridge=%d trunc=%d; H5 cE=%.4e erg/g "
          "(std/kT %.3f, n=%d); overlap p-mismatch max %.2e"
          % (f_melt, fit_rho[1], st["cols"], st["zone1"], st["blend"],
             st["bridge"], st["blend_trunc"], st["align_cE"],
             st["align_std_kT"], st["align_n"], st["p_mismatch_max"]))
    st["rho_max"] = fit_rho[1]

    p0 = p
    p2, e2, _mask, mx2 = condition_surface(10.0 ** np.asarray(lrho, float),
                                           10.0 ** np.asarray(lT, float),
                                           p, e)
    sig = (np.abs(p2 - p0) > 1e-6 * np.maximum(np.abs(p0), 1e-300)) | (p0 <= 0.0)
    n_new = int((sig & ~band).sum())
    band = band | sig
    hull = np.where(sig, 0.0, hull)
    print("crossover-2 (post-blend): %d maxwell / %d flat / %d ramp "
          "isotherms; %d newly banded cells"
          % (mx2["n_maxwell"], mx2["n_flat"], mx2["n_ramp"], n_new))
    # NOTE (measured, 2963): a post-crossover per-column re-grade of the
    # gap columns' T-structure was tried and removed — monotonise_rho's
    # cummax against the vapor-side fill overwrites any per-column
    # T-grading at rho < rho0, so the emitted surface is identical either
    # way (2-D monotonicity squeeze). A genuinely T-graded gap needs the
    # deferred R1 vapor-branch model, not a patch here.
    st.update(f_melt=f_melt, fit=dict(vinet=vinet, rms=rms, n=n_fit),
              model=rep, cx2_new_band=n_new, cx2=mx2)
    return p2, e2, hull, band, st
