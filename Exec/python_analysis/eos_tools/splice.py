"""Tier-1 splice: assemble the wide-range deuterium table (plan SS3).

Sources (cold -> hot), weights, and seams (v2, doc/eos_creation_plan.md):

    cold model | H-REOS.3(D) | FPEOS(D) | ideal plasma
    seam 1: lT = 3.2 +- 0.2      (cold <-> REOS.3)
    seam 2: lT = 5.5 +- 0.25     (REOS.3 <-> FPEOS)
    seam 3: lT = 7.6 +- 0.2      (FPEOS <-> ideal plasma)
    rho hand-off: cold model only below ~1 g/cc (lrho ramp 0.0 +- 0.1) —
    beyond the LLNL cold-curve fit range REOS.3's measured 60 K+ data wins.

Nominal weight partition (sums to 1 identically):

    A_cold=(1-s1); A_reos=s1(1-s2); A_fp=s2(1-s3); A_ip=s3
    W = [A_cold*rc, A_reos + A_cold*(1-rc), A_fp, A_ip],  rc = cold-capable

Per-cell hull handling: weights of sources whose hull is 0 (or that were
not evaluated there) are dropped and the rest renormalised; a cell keeping
< 95% of its nominal weight is marked hull = 0 (its value comes from the
best available fallback, e.g. the saturated-liquid fill inside the dome or
the REOS.3 60 K row below 60 K at high rho).

Energy alignment (splice-plan §5.1): constant e-offsets chained cold ->
REOS.3 -> FPEOS -> ideal plasma, measured in each seam band over cells
where both partners are in-hull; offset constancy (std / local kT-scale)
is a QA gate. Pressures are never shifted.

Tier 1 blends p and e directly (T-only seams; the omitted Helmholtz
correction is measured by the Maxwell-residual QA gate, plan §1.2).
"""

import gzip
import os
import shutil

import numpy as np

from .condition import condition, inverse_maps, shift_energy
from .constants import KB, M_D
from .formats.fpeos import fpeos_to_deuterium, read_fpeos
from .formats.reos3 import read_reos3, reos3_to_deuterium
from .grids import regrid_onto
from .models.composite import LogRamp
from . import sources as _sources

# --- deuterium splice constants -------------------------------------------
GRID = dict(lrho=(-4.0, 3.0, 384), lT=(1.25, 9.0, 384))
# Seam 1 sits BELOW the cold model's internal CoolProp->rotor-gas ramp
# (505-600 K): the C++ self-test on the first emitted table caught a 24%
# p(T) dip at (0.2 g/cc, ~533 K) — the ideal rotor gas blending against a
# dense supercritical fluid. The rotor piece was an iFPEOS-floor (800 K)
# crutch; H-REOS.3 reaches 60 K, so the table hands to it at ~355 K and the
# rotor gas never carries splice weight. Band top (+1.5 hw, the C1 metric
# edge) stays under CoolProp's 600 K ceiling.
S1 = LogRamp(2.55, 0.15)
S2 = LogRamp(5.5, 0.25)
S3 = LogRamp(7.6, 0.20)
RHOC_C = -0.62                     # rho-handoff center (log10; ~0.24 g/cc)
RHOC_HW = 0.12                     # rho-handoff ramp halfwidth (dex)


def rho_capable(R, TT):
    """Cold-model density weight rc(rho): constant hand-off at ~0.24 g/cc.

    Two measured iterations fixed this value (see eos_creation_plan.md):
    - hand-off at 1 g/cc exposed the cold model at (0.3-0.5 g/cc, 1-2.5 kK)
      — QEOS solid far above melt vs REOS.3's dense dissociating fluid —
      and the local e-mismatch put a 0.43 compression discontinuity on the
      principal Hugoniot at ~10 GPa (right where dissociation steepens it);
    - a T-dependent (sliding) hand-off over-corrected: |dc/dlT| ~ 2.7 made
      a diagonal weight transition far sharper in T than the seam ramps
      (seam-1 slope-jump 0.09 -> 0.82).
    Constant hand-off at 0.24 g/cc keeps the reference state (0.171 g/cc)
    cold-owned with margin and cedes the mismatch region to REOS.3 at all
    T. Cost: the compressed-cryo corner (0.24-1 g/cc below 60 K) drops to
    hull-0 REOS-60 K fill — off every shock path, e error ~ cv*(60-20 K)
    ~ 1e-5 relative.
    """
    x = np.log10(np.asarray(R, float))
    K = 1.47
    return (1.0 - 0.5 * (1.0 + np.tanh(K * (x - RHOC_C) / RHOC_HW))) \
        * np.ones_like(np.asarray(TT, float))
P_FLOOR = 1.0e3                     # barye (tension clip, decision 4)
W_EPS = 1e-4                        # nominal weights below this are zeroed
HULL_KEEP = 0.95                    # min kept nominal weight for hull = 1
SRC_NAMES = ("cold_model_d2", "H-REOS.3(D)", "FPEOS(D)", "ideal_plasma")

REF_STATE = dict(rho0=0.171, T0=20.0)  # cryogenic liquid D2


def raw_path(*parts):
    return os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                        "data", "raw", *parts)


def nominal_weights(R, TT):
    """(4, nr, nt) nominal weight stack; sums to 1."""
    s1, s2, s3 = S1.w(TT), S2.w(TT), S3.w(TT)
    rc = rho_capable(R, TT)
    A_cold, A_reos = (1 - s1), s1 * (1 - s2)
    A_fp, A_ip = s2 * (1 - s3), s3
    W = np.stack([A_cold * rc, A_reos + A_cold * (1 - rc), A_fp, A_ip])
    W[W < W_EPS] = 0.0
    return W / W.sum(axis=0)


def _subrect(mask):
    """Bounding index slices of a 2-D mask (assumes non-empty)."""
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    return slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1)


def build_sources(lrho, lT, W, use_coolprop=True, fpeos_src=None):
    """Evaluate/regrid all four sources onto the grid.

    Returns dict name -> {p, e, hull, avail} (full-shape arrays)."""
    nr, nt = len(lrho), len(lT)
    R = 10.0 ** lrho[:, None] * np.ones((1, nt))
    TT = 10.0 ** lT[None, :] * np.ones((nr, 1))
    out = {}

    # cold composite (analytic; evaluate only where it carries weight)
    from .materials.deuterium import DeuteriumColdModel
    cm = DeuteriumColdModel(use_coolprop=use_coolprop)
    m = W[0] > 0.0
    sr, sc = _subrect(m)
    p = np.zeros((nr, nt))
    e = np.zeros((nr, nt))
    hull = np.zeros((nr, nt))
    avail = np.zeros((nr, nt), bool)
    cd = cm.eval(R[sr, sc], TT[sr, sc])
    p[sr, sc], e[sr, sc] = cd["p"], cd["e"]
    hull[sr, sc] = cd["hull"]
    avail[sr, sc] = True
    out["cold"] = dict(p=p, e=e, hull=hull, avail=avail, model=cm)

    # H-REOS.3 (D basis)
    reos = reos3_to_deuterium(read_reos3(raw_path("reos3",
                                                  "table2_HREOS3.dat")))
    p, e, hull, loo = regrid_onto(reos, lrho, lT)
    out["reos"] = dict(p=p, e=e, hull=hull,
                       avail=np.ones((nr, nt), bool), loo=loo)

    # FPEOS (D basis) — read from the committed tarball unless given a path
    if fpeos_src is None:
        import tarfile
        import tempfile
        with tarfile.open(raw_path("fpeos_10-26-25.tar.gz"), "r:gz") as tf:
            txt = tf.extractfile("FPEOS/H_EOS_09-18-20.txt").read().decode()
        with tempfile.NamedTemporaryFile("w", suffix=".txt",
                                         delete=False) as f:
            f.write(txt)
            fpeos_src = f.name
        pts = fpeos_to_deuterium(read_fpeos(fpeos_src))
        os.unlink(fpeos_src)
    else:
        pts = fpeos_to_deuterium(read_fpeos(fpeos_src))
    p, e, hull, loo = regrid_onto(pts, lrho, lT)
    out["fpeos"] = dict(p=p, e=e, hull=hull,
                        avail=np.ones((nr, nt), bool), loo=loo)

    # ideal plasma (fast tabulated evaluator). Evaluated well below its
    # weight region (lT >= 5) so the fallback ladder can use it wherever
    # FPEOS has no hull (rho below the FPEOS floor stays fully ionized and
    # weakly coupled at these temperatures).
    from .models.ideal_plasma import IdealPlasmaFast
    ip = IdealPlasmaFast(M_D, Z=1.0)
    m = (W[3] > 0.0) | (lT[None, :] >= 5.0)
    p = np.zeros((nr, nt))
    e = np.zeros((nr, nt))
    hull = np.zeros((nr, nt))
    avail = np.zeros((nr, nt), bool)
    if m.any():
        sr, sc = _subrect(m * np.ones((nr, 1), bool))
        p[sr, sc] = ip.p(R[sr, sc], TT[sr, sc])
        e[sr, sc] = ip.e(R[sr, sc], TT[sr, sc])
        hull[sr, sc] = 1.0
        avail[sr, sc] = True
    out["ip"] = dict(p=p, e=e, hull=hull, avail=avail)
    return out


def align_chain(srcs, W, lrho, lT):
    """Chained constant e-offsets, applied in place. Returns stats list."""
    TT = 10.0 ** lT[None, :] * np.ones((len(lrho), 1))
    kT = 1.5 * KB * TT / M_D  # local energy scale, erg/g
    order = ["cold", "reos", "fpeos", "ip"]
    bands = [np.abs(lT[None, :] - S1.c) <= S1.d,
             np.abs(lT[None, :] - S2.c) <= S2.d,
             np.abs(lT[None, :] - S3.c) <= S3.d]
    stats = []
    for k, (lo, hi) in enumerate(zip(order[:-1], order[1:])):
        a, b = srcs[lo], srcs[hi]
        sel = (bands[k] * np.ones_like(TT, bool)) \
            & (a["hull"] > 0.5) & (b["hull"] > 0.5) \
            & a["avail"] & b["avail"]
        if k == 0:  # seam 1 only where the cold model actually holds
            R1 = 10.0 ** lrho[:, None] * np.ones((1, len(lT)))
            sel &= rho_capable(R1, TT) > 0.5
        de = a["e"][sel] - b["e"][sel]
        cE = float(np.mean(de))
        b["e"] = b["e"] + cE
        stats.append(dict(pair="%s->%s" % (hi, lo), n=int(sel.sum()), cE=cE,
                          std_over_kT=float(np.mean(np.abs(de - cE)
                                                    / kT[sel]))))
    return stats


def blend_all(srcs, W, lT):
    """Hull-aware weighted blend. Returns (p, e, hull, kept_weight).

    Cells where every weighted source is invalid fall back to a single
    source by a physics-priority ladder — never to regrid FILL values of
    a data source outside its hull (the first emitted table used the
    FPEOS nearest-fill below its density floor, which planted a large
    non-monotone p/e bump across every dilute row; the C++ self-test
    round-trips caught it):
        lT >= 6:  ideal plasma  ->  REOS.3 fill  ->  cold model
        lT <  6:  cold model    ->  REOS.3 fill
    (REOS.3 fill remains acceptable: its only fallback uses are the
    below-60 K/high-rho corner and the dome edges, where the nearest row
    is physically adjacent.)
    """
    names = ["cold", "reos", "fpeos", "ip"]
    valid = np.stack([(srcs[n]["hull"] > 0.5) & srcs[n]["avail"]
                      for n in names])
    P = np.stack([srcs[n]["p"] for n in names])
    E = np.stack([srcs[n]["e"] for n in names])

    Wv = W * valid
    kept = Wv.sum(axis=0)
    hull = kept >= HULL_KEEP
    Wn = Wv / np.maximum(kept, 1e-300)
    p = (Wn * P).sum(axis=0)
    e = (Wn * E).sum(axis=0)

    # fallback ladder for kept == 0
    need = kept <= 0.0
    if need.any():
        hot = need & (lT[None, :] >= 6.0)
        cold_av = srcs["cold"]["avail"]
        ip_av = srcs["ip"]["avail"]
        # priority stacks (first available wins)
        for region, order in ((hot, ["ip", "reos", "cold"]),
                              (need & ~hot, ["cold", "reos"])):
            left = region.copy()
            for nsrc in order:
                av = srcs[nsrc]["avail"] if nsrc != "reos" \
                    else np.ones_like(cold_av)
                take = left & av
                p = np.where(take, srcs[nsrc]["p"], p)
                e = np.where(take, srcs[nsrc]["e"], e)
                left &= ~take
    return p, e, hull.astype(float), kept


def monotonise_T(F, rel_eps=1e-12):
    """Enforce strictly nondecreasing F(T) along every rho-row.

    The v1 inversion contract (single-branch invert_T_from_e/p) needs
    monotone value surfaces; the conditioning already floors cv > 0, so
    this makes the values consistent with the derivative blocks. Measured
    magnitudes are small (worst p dip 1.2e-2 relative, e 3.6e-4, in the
    REOS.3 near-melt rows). Returns (F', n_changed, max_rel_change).
    """
    Fm = np.maximum.accumulate(F, axis=1)
    # strictify flats with a monotone epsilon ramp (invisible at 1e-12;
    # sign-safe: added, not multiplied — e is negative before e_shift)
    scale = np.maximum(np.abs(Fm), 1e-300)
    Fm = np.maximum.accumulate(
        Fm + rel_eps * scale * np.arange(F.shape[1])[None, :], axis=1)
    tol = 10.0 * rel_eps * F.shape[1] * np.maximum(np.abs(F), 1e-300)
    changed = (Fm - F) > tol
    max_rel = float(((Fm - F) / np.maximum(np.abs(F), 1e-300)).max())
    return Fm, int(changed.sum()), max_rel


def splice_deuterium(n_rho=None, n_T=None, use_coolprop=True,
                     p_floor=P_FLOOR):
    """Run the full pipeline. Returns dict with everything QA needs."""
    lr0, lr1, nr = GRID["lrho"]
    lt0, lt1, nt = GRID["lT"]
    nr, nt = n_rho or nr, n_T or nt
    lrho = np.linspace(lr0, lr1, nr)
    lT = np.linspace(lt0, lt1, nt)
    R = 10.0 ** lrho[:, None] * np.ones((1, nt))
    TT = 10.0 ** lT[None, :] * np.ones((nr, 1))

    W = nominal_weights(R, TT)
    srcs = build_sources(lrho, lT, W, use_coolprop=use_coolprop)
    align_stats = align_chain(srcs, W, lrho, lT)
    p, e, hull, kept = blend_all(srcs, W, lT)

    # tension clip (decision 4): keeps the C++ in-hull p_min > 0 contract
    clip = p < p_floor
    n_clip = int(clip.sum())
    p = np.where(clip, p_floor, p)
    hull = np.where(clip, 0.0, hull)

    # monotone-in-T value surfaces (single-branch inversion contract)
    p, n_mono_p, mrel_p = monotonise_T(p)
    e, n_mono_e, mrel_e = monotonise_T(e)
    mono_stats = dict(n_p=n_mono_p, rel_p=mrel_p, n_e=n_mono_e, rel_e=mrel_e)

    e, e_shift = shift_energy(e, hull)
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, p, e)
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks = dict(p=p, e=e, dpdT=dpdT, dpdrho=dpdrho, cv=cv, dedrho=dedrho,
                  T_of_e=T_of_e, T_of_p=T_of_p, hull=hull)
    return dict(lrho=lrho, lT=lT, blocks=blocks, inv_axes={"le": le,
                "lp": lp}, e_shift=e_shift, cond_stats=stats,
                align_stats=align_stats, srcs=srcs, W=W, kept=kept,
                n_clip=n_clip, p_floor=p_floor, mono_stats=mono_stats)


def write_spliced(res, out_path, generator_line):
    """Emit .eostab + .gz with splice provenance."""
    from .formats.eostab import write_eostab
    st = res["cond_stats"]
    prov = [
        ("material", "D"),
        ("source", "spliced: cold_model_d2(CoolProp/Richardson2014 + "
         "LLNL-686936 Vinet + Slater-Debye) | H-REOS.3(D-scaled, ApJS 215 "
         "21) | FPEOS(D-scaled, PRE 103 013203) | ideal-plasma; seams "
         "lT=%.2f/%.2f/%.2f; see doc/eos_creation_plan.md + "
         "data/raw/sources.yaml" % (S1.c, S2.c, S3.c)),
        ("generator", generator_line),
        ("composition", "A=2.014 Z=1"),
        ("units", "cgs"),
        ("e_shift", "%.10e" % res["e_shift"]),
        ("conditioning", "cv_floor=%.4e cv_floored=%d monotonised=%d "
         "maxwell=%s tension_clip=%d monoT_p=%d(%.1e) monoT_e=%d(%.1e) "
         "align_e=%s"
         % (st["cv_floor"], st["cv_floored"], st["monotonised"],
            st["maxwell"], res["n_clip"],
            res["mono_stats"]["n_p"], res["mono_stats"]["rel_p"],
            res["mono_stats"]["n_e"], res["mono_stats"]["rel_e"],
            ",".join("%.4e" % a["cE"] for a in res["align_stats"]))),
    ]
    write_eostab(out_path, prov, res["lrho"], res["lT"], res["blocks"],
                 res["inv_axes"])
    with open(out_path, "rb") as fin, gzip.open(out_path + ".gz", "wb",
                                                compresslevel=9) as fout:
        shutil.copyfileobj(fin, fout)
    print("wrote %s.gz" % out_path)
