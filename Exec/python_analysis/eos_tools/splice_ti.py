"""Titanium splice (Ti SS3-prime): solid model | Saha | FD ideal plasma.

Three sources on the D-splice pattern (T-only tanh seams, hull-aware
weights, constant-cE energy alignment, monotone-in-T conditioning):

    solid  = TitaniumColdModel (Vinet + Slater-Debye + Sommerfeld)
    saha   = SahaModel with the NIST Ti I-XXII energies (ideal mixture;
             also the low-density vapor limit down to the T floor)
    ip     = IdealPlasmaFast(Z = 22) — the degeneracy-correct hot cap

    seam 1: lT = 4.50 +- 0.25  (solid <-> Saha)
    seam 3: lT = 8.40 +- 0.20  (Saha <-> ideal plasma; Zbar ~ 22 there)
    rho hand-off: solid only above ~3.2 g/cc (log ramp 0.505 +- 0.12)

**Ti v0 accuracy statement** (doc/eos_splice_plan.md §4 was explicit that
this band is genuinely new work): the seam-1 band and the compressed WDM
region carry uncontrolled error — the solid model extrapolates Debye +
Sommerfeld beyond their comfort and the Saha side has no cohesion,
degeneracy or continuum lowering. The in-band alignment scatter and C1
metrics are therefore REPORTS for Ti v0, not hard gates; hard gates are
structural (convexity, monotone surfaces, finiteness, Hugoniot smoothness,
C++ round-trips). Filling the band needs OFMD/ML-MD data (manifest:
'ti-mlmd-melt', pending manual).

Known-physical hull-0 regions: the expanded-condensed ("spall") wedge
(0.02-3.2 g/cc below 4 kK — ideal-vapor values, two-phase in reality) and
the solid tension clip below the ambient density at low T.
"""

import numpy as np

from .condition import condition, inverse_maps, shift_energy
from .constants import KB
from .models.composite import LogRamp
from .splice import (HULL_KEEP, W_EPS, monotonise_T)

GRID_TI = dict(lrho=(-3.0, 1.7, 384), lT=(1.25, 9.0, 384))
S1_TI = LogRamp(4.50, 0.25)
S3_TI = LogRamp(8.40, 0.20)
RHOC_TI = LogRamp(np.log10(3.2), 0.12)
P_FLOOR_TI = 1.0e3  # barye

# expanded-condensed wedge: flagged hull-0 (two-phase in reality)
WEDGE_RHO = (0.02, 3.2)   # g/cc
WEDGE_TMAX = 4.0e3        # K

REF_STATE_TI = dict(rho0=4.51, T0=293.0)


def nominal_weights_ti(R, TT):
    s1, s3 = S1_TI.w(TT), S3_TI.w(TT)
    # the solid is capable at HIGH rho: its density weight is the ramp
    # itself (1 above ~3.2 g/cc), unlike deuterium's low-rho cold model
    w_solid_rho = RHOC_TI.w(R)
    A_solid, A_saha, A_ip = (1 - s1), s1 * (1 - s3), s3
    W = np.stack([A_solid * w_solid_rho,
                  A_saha + A_solid * (1.0 - w_solid_rho),
                  A_ip])
    W[W < W_EPS] = 0.0
    return W / W.sum(axis=0)


def _subrect(mask):
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    return slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1)


def build_sources_ti(lrho, lT, W):
    from .materials.titanium import M_TI, TitaniumColdModel
    from .models.ideal_plasma import IdealPlasmaFast
    from .models.saha import titanium_saha

    nr, nt = len(lrho), len(lT)
    R = 10.0 ** lrho[:, None] * np.ones((1, nt))
    TT = 10.0 ** lT[None, :] * np.ones((nr, 1))
    out = {}

    solid = TitaniumColdModel()
    m = W[0] > 0.0
    p = np.zeros((nr, nt))
    e = np.zeros((nr, nt))
    avail = np.zeros((nr, nt), bool)
    if m.any():
        sr, sc = _subrect(m)
        sd = solid.eval(R[sr, sc], TT[sr, sc])
        p[sr, sc], e[sr, sc] = sd["p"], sd["e"]
        avail[sr, sc] = True
    out["solid"] = dict(p=p, e=e, hull=avail.astype(float), avail=avail,
                        model=solid)

    saha = titanium_saha()
    sd = saha.eval(R, TT)
    out["saha"] = dict(p=sd["p"], e=sd["e"], zbar=sd["zbar"],
                       hull=np.ones((nr, nt)),
                       avail=np.ones((nr, nt), bool), model=saha)

    ip = IdealPlasmaFast(M_TI, Z=22.0)
    m = (W[2] > 0.0) | (lT[None, :] >= 7.4)
    p = np.zeros((nr, nt))
    e = np.zeros((nr, nt))
    avail = np.zeros((nr, nt), bool)
    sr, sc = _subrect(m * np.ones((nr, 1), bool))
    p[sr, sc] = ip.p(R[sr, sc], TT[sr, sc])
    e[sr, sc] = ip.e(R[sr, sc], TT[sr, sc])
    avail[sr, sc] = True
    out["ip"] = dict(p=p, e=e, hull=avail.astype(float), avail=avail)
    return out


def splice_titanium(n_rho=None, n_T=None, p_floor=P_FLOOR_TI):
    lr0, lr1, nr = GRID_TI["lrho"]
    lt0, lt1, nt = GRID_TI["lT"]
    nr, nt = n_rho or nr, n_T or nt
    lrho = np.linspace(lr0, lr1, nr)
    lT = np.linspace(lt0, lt1, nt)
    R = 10.0 ** lrho[:, None] * np.ones((1, nt))
    TT = 10.0 ** lT[None, :] * np.ones((nr, 1))

    W = nominal_weights_ti(R, TT)
    srcs = build_sources_ti(lrho, lT, W)

    # energy alignment chain: saha -> solid (seam-1 band, solid side),
    # then ip -> aligned saha (seam-3 band)
    kT = 1.5 * KB * TT / srcs["saha"]["model"].m
    align_stats = []
    band1 = (np.abs(lT[None, :] - S1_TI.c) <= S1_TI.d) \
        & srcs["solid"]["avail"] & (RHOC_TI.w(R) > 0.5)
    de = srcs["solid"]["e"][band1] - srcs["saha"]["e"][band1]
    cE1 = float(np.mean(de))
    srcs["saha"]["e"] = srcs["saha"]["e"] + cE1
    align_stats.append(dict(pair="saha->solid", n=int(band1.sum()), cE=cE1,
                            std_over_kT=float(np.mean(np.abs(de - cE1)
                                                      / kT[band1]))))
    band3 = (np.abs(lT[None, :] - S3_TI.c) <= S3_TI.d) & srcs["ip"]["avail"]
    de = srcs["saha"]["e"][band3] - srcs["ip"]["e"][band3]
    cE3 = float(np.mean(de))
    srcs["ip"]["e"] = srcs["ip"]["e"] + cE3
    align_stats.append(dict(pair="ip->saha", n=int(band3.sum()), cE=cE3,
                            std_over_kT=float(np.mean(np.abs(de - cE3)
                                                      / kT[band3]))))

    # hull-aware blend + physics-priority fallback (saha always available)
    names = ["solid", "saha", "ip"]
    valid = np.stack([(srcs[n]["hull"] > 0.5) & srcs[n]["avail"]
                      for n in names])
    P = np.stack([srcs[n]["p"] for n in names])
    E = np.stack([srcs[n]["e"] for n in names])
    Wv = W * valid
    kept = Wv.sum(axis=0)
    hull = (kept >= HULL_KEEP).astype(float)
    Wn = Wv / np.maximum(kept, 1e-300)
    p = (Wn * P).sum(axis=0)
    e = (Wn * E).sum(axis=0)
    need = kept <= 0.0
    if need.any():  # only possible where solid+ip unavailable -> saha
        p = np.where(need, srcs["saha"]["p"], p)
        e = np.where(need, srcs["saha"]["e"], e)
        hull = np.where(need, 0.0, hull)

    # expanded-condensed wedge: values stand, trustworthiness flagged off
    wedge = (R > WEDGE_RHO[0]) & (R < WEDGE_RHO[1]) & (TT < WEDGE_TMAX)
    hull = np.where(wedge, 0.0, hull)

    clip = p < p_floor
    n_clip = int(clip.sum())
    p = np.where(clip, p_floor, p)
    hull = np.where(clip, 0.0, hull)

    p, n_mono_p, mrel_p = monotonise_T(p)
    e, n_mono_e, mrel_e = monotonise_T(e)

    e, e_shift = shift_energy(e, hull)
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, p, e)
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks = dict(p=p, e=e, dpdT=dpdT, dpdrho=dpdrho, cv=cv, dedrho=dedrho,
                  T_of_e=T_of_e, T_of_p=T_of_p, hull=hull)
    return dict(lrho=lrho, lT=lT, blocks=blocks,
                inv_axes={"le": le, "lp": lp}, e_shift=e_shift,
                cond_stats=stats, align_stats=align_stats, srcs=srcs, W=W,
                kept=kept, n_clip=n_clip, p_floor=p_floor,
                mono_stats=dict(n_p=n_mono_p, rel_p=mrel_p, n_e=n_mono_e,
                                rel_e=mrel_e))


def write_spliced_ti(res, out_path, generator_line):
    import gzip
    import shutil

    from .formats.eostab import write_eostab
    st = res["cond_stats"]
    prov = [
        ("material", "Ti"),
        ("source", "spliced Ti v0: solid(LLNL-686936 Vinet + Slater-Debye "
         "+ Sommerfeld) | Saha(NIST Ti I-XXII, ideal) | FD ideal plasma "
         "Z=22; seams lT=%.2f/%.2f, rho hand-off 3.2 g/cc. WDM band "
         "accuracy UNCONTROLLED pending OFMD/ML-MD data — see "
         "doc/eos_creation_plan.md + data/raw/sources.yaml"
         % (S1_TI.c, S3_TI.c)),
        ("generator", generator_line),
        ("composition", "A=47.867 Z=22"),
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
