#!/usr/bin/env python3
"""Offline EOS table conditioning tool for Cerberus tabulated-EOS (plan W3).

Produces canonical ``.eostab`` files from raw source data (FPEOS first),
plus QA plots. All curation happens HERE, offline and inspectable; the C++
reader (Stage 2, MFP_eos_table) parses only the canonical format below.

.eostab format spec v1 — FROZEN 2026-07-07
(mirror of Exec/testing/EOS-Table/README.md; do not change one without the
other):

    EOSTAB 1                          # magic + format version
    # provenance: free-form 'key: value' lines, '#' comments allowed
    material:    D
    source:      <origin, retrieval date>
    generator:   eos_table_prep.py <git-sha>, run <date>
    composition: A=2.014 Z=1
    units:       cgs                  # values stored DIMENSIONAL (CGS):
                                      #   rho g/cc, T K, p erg/cc, e erg/g
                                      # Cerberus nondimensionalises at load
    e_shift:     <erg/g>              # constant added to raw specific
                                      # internal energy so min(e) > 0 on the
                                      # hull (zero-point is arbitrary; must
                                      # be used consistently, and is)
    conditioning: cv_floor=<v> monotonised=<n> cv_floored=<n> maxwell=<status>
    # grid (uniform in log10); axes reconstructed by reader, never stored
    grid: n_rho=<Nr> n_T=<Nt>
    lrho: <log10 rho_min> <log10 rho_max>
    lT:   <log10 T_min>   <log10 T_max>
    # inverse-map axes (required iff T_of_e / T_of_p blocks present)
    le:   <log10 e_min> <log10 e_max>  n_e=<Ne>
    lp:   <log10 p_min> <log10 p_max>  n_p=<Np>
    # data blocks: 'block: <name>' then the values, whitespace separated.
    # storage order: i_rho outer, j inner  (idx = i*n_T + j; inverse maps
    # idx = i*n_e + j / i*n_p + j)
    block: p          # pressure
    block: e          # specific internal energy (shifted)
    block: dpdT       # (dP/dT)|rho     smoothed  -- outputs only
    block: dpdrho     # (dP/drho)|T     smoothed/monotonised
    block: cv        # (de/dT)|rho     smoothed, floored > 0
    block: dedrho     # (de/drho)|T     smoothed
    block: T_of_e     # T(rho, e) inverse map  -- Newton seeds (plan D5/D8)
    block: T_of_p     # T(rho, p) inverse map
    block: hull       # REQUIRED: 1.0 inside source hull, 0.0 filled cell

Commands:
    synthetic  — analytic ideal-gas table (Stage-2 tier-1 gate + fallback)
    fpeos      — ingest FPEOS H table, isotope-scale to D, regrid,
                 condition, write, QA
    qa         — QA plots for an existing .eostab

Examples (from Exec/testing/EOS-Table/):
    python3 ../../python_analysis/eos_table_prep.py synthetic \
        --out data/ideal_synthetic.eostab
    python3 ../../python_analysis/eos_table_prep.py fpeos \
        --src data/raw/FPEOS/H_EOS_09-18-20.txt \
        --out data/D_fpeos.eostab --qa qa
"""

import argparse
import datetime
import os
import subprocess

import numpy as np
from scipy.interpolate import PchipInterpolator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --- physical constants (CGS) ------------------------------------------------
HA_ERG = 4.3597447222071e-11  # Hartree -> erg
AMU_G = 1.66053906660e-24  # atomic mass unit -> g
M_H = 1.00782503207 * AMU_G  # hydrogen-atom mass
M_D = 2.01355321271 * AMU_G  # deuterium-atom mass
KB = 1.380649e-16  # Boltzmann, erg/K
GPA_CGS = 1.0e10  # GPa -> erg/cc (dyn/cm^2)


def git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
        ).decode().strip()
    except Exception:
        return "unknown"


# =============================================================================
# writer
# =============================================================================


def write_eostab(path, prov, lrho, lT, blocks, inv_axes):
    """Write a .eostab per the frozen spec (see module docstring).

    prov: list of (key, value) provenance lines, in order.
    lrho/lT: 1-D log10 axes (uniform).
    blocks: dict name -> 2-D array; forward blocks shaped (Nr, Nt),
            T_of_e (Nr, Ne), T_of_p (Nr, Np).
    inv_axes: dict with 'le' (1-D log10 e axis), 'lp' (1-D log10 p axis).
    """
    nr, nt = len(lrho), len(lT)
    order = ["p", "e", "dpdT", "dpdrho", "cv", "dedrho", "T_of_e", "T_of_p", "hull"]
    with open(path, "w") as f:
        f.write("EOSTAB 1\n")
        for k, v in prov:
            f.write("%s: %s\n" % (k, v))
        f.write("grid: n_rho=%d n_T=%d\n" % (nr, nt))
        f.write("lrho: %.12g %.12g\n" % (lrho[0], lrho[-1]))
        f.write("lT: %.12g %.12g\n" % (lT[0], lT[-1]))
        le, lp = inv_axes["le"], inv_axes["lp"]
        f.write("le: %.12g %.12g n_e=%d\n" % (le[0], le[-1], len(le)))
        f.write("lp: %.12g %.12g n_p=%d\n" % (lp[0], lp[-1], len(lp)))
        for name in order:
            arr = np.asarray(blocks[name])
            assert arr.shape[0] == nr, (name, arr.shape)
            assert np.all(np.isfinite(arr)), "non-finite values in block " + name
            f.write("block: %s\n" % name)
            flat = arr.ravel(order="C")  # idx = i*ncol + j
            for i in range(0, flat.size, 6):
                f.write(" ".join("%.10e" % v for v in flat[i:i + 6]) + "\n")
    print("wrote %s  (%d x %d, %d blocks)" % (path, nr, nt, len(order)))


def read_eostab(path):
    """Minimal reader (QA use; prototype of the C++ parser)."""
    prov, blocks = {}, {}
    with open(path) as f:
        tok = f.read().split("\n")
    assert tok[0].startswith("EOSTAB 1"), "bad magic"
    i, cur = 1, None
    vals = []
    meta = {}
    while i < len(tok):
        line = tok[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if line.startswith("block:"):
            if cur:
                blocks[cur] = np.array(vals)
            cur, vals = line.split()[1], []
        elif cur is not None:
            vals.extend(float(x) for x in line.split())
        else:
            k, v = line.split(":", 1)
            prov[k.strip()] = v.strip()
    if cur:
        blocks[cur] = np.array(vals)
    g = dict(p.split("=") for p in prov["grid"].split())
    nr, nt = int(g["n_rho"]), int(g["n_T"])
    a, b = (float(x) for x in prov["lrho"].split())
    meta["lrho"] = np.linspace(a, b, nr)
    a, b = (float(x) for x in prov["lT"].split())
    meta["lT"] = np.linspace(a, b, nt)
    for ax, nkey in (("le", "n_e"), ("lp", "n_p")):
        parts = prov[ax].split()
        n = int(dict(p.split("=") for p in parts if "=" in p)[nkey])
        meta[ax] = np.linspace(float(parts[0]), float(parts[1]), n)
    for name in blocks:
        ncol = {"T_of_e": len(meta["le"]), "T_of_p": len(meta["lp"])}.get(name, nt)
        blocks[name] = blocks[name].reshape(nr, ncol)
    return prov, meta, blocks


# =============================================================================
# FPEOS ingest
# =============================================================================


def read_fpeos(path):
    """Parse an FPEOS *_EOS_*.txt element table.

    Line format:
    f= H N= 1 rho[g/cc]= v V[A^3]= v T[K]= v P[GPa]= v err E[Ha]= v err
    Returns dict of 1-D arrays (per-atom energies, hydrogen-mass densities)
    and the Hugoniot initial condition from the header.
    """
    rho, T, P, E = [], [], [], []
    e0_ha, v0_a3 = None, None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "E0[Ha]=" in line:
                    e0_ha = float(line.split("E0[Ha]=")[1].split()[0])
                if "V0[A^3]=" in line:
                    v0_a3 = float(line.split("V0[A^3]=")[1].split()[0])
                continue
            t = line.split()
            if len(t) < 16 or t[0] != "f=":
                continue
            rho.append(float(t[5]))
            T.append(float(t[9]))
            P.append(float(t[11]))
            E.append(float(t[14]))
    return {
        "rho": np.array(rho), "T": np.array(T),
        "P_GPa": np.array(P), "E_Ha": np.array(E),
        "E0_Ha": e0_ha, "V0_A3": v0_a3, "n_atoms_per_formula": 1,
    }


def fpeos_to_deuterium(raw):
    """Isotope-scale the per-atom hydrogen table to deuterium (CGS).

    Equal nuclear number density and T: rho_D = rho_H * m_D/m_H,
    P unchanged, e[erg/g] = E[Ha/atom]*HA_ERG/m_D.
    """
    s = M_D / M_H
    pts = {
        "rho": raw["rho"] * s,
        "T": raw["T"].copy(),
        "p": raw["P_GPa"] * GPA_CGS,
        "e": raw["E_Ha"] * HA_ERG / M_D,
    }
    # Hugoniot initial condition, same scaling (V0 is per atom)
    pts["hug_rho0"] = M_D / (raw["V0_A3"] * 1e-24)
    pts["hug_e0"] = raw["E0_Ha"] * HA_ERG / M_D
    pts["hug_p0"] = 0.0
    return pts


# =============================================================================
# regrid + condition
# =============================================================================


def isochores(pts):
    """Group scattered points into per-density isochores, T-sorted."""
    out = []
    for r in np.unique(pts["rho"]):
        m = pts["rho"] == r
        o = np.argsort(pts["T"][m])
        out.append((r, pts["T"][m][o], pts["p"][m][o], pts["e"][m][o]))
    return out


def regrid(pts, n_rho, n_T):
    """Hull-aware two-pass PCHIP regrid onto a log-uniform (rho, T) grid.

    Pass a: along each source isochore in log T (within its own T range).
    Pass b: along log rho at each target T (within the covered rho range).
    Cells never covered are filled with the nearest hull value and marked
    0 in the hull mask. Refuses to extrapolate anywhere.
    Returns (lrho, lT, p, e, hull, loo) with p,e shaped (n_rho, n_T);
    loo = leave-one-out regrid-error stats (review amendment: measure,
    don't assume, on a coarse source).
    """
    iso = isochores(pts)
    lrho = np.linspace(np.log10(min(r for r, *_ in iso)),
                       np.log10(max(r for r, *_ in iso)), n_rho)
    lT = np.linspace(np.log10(pts["T"].min()), np.log10(pts["T"].max()), n_T)

    # pass a: isochore -> target T nodes (log-log for p; log T, linear e)
    partial = []  # (log10 rho, mask, p_col, e_col)
    loo_err = []
    for r, Ts, Ps, Es in iso:
        if len(Ts) < 4:
            continue  # too sparse to shape a curve; drop (recorded via hull)
        lt = np.log10(Ts)
        fp = PchipInterpolator(lt, np.log10(Ps))
        fe = PchipInterpolator(lt, Es)
        m = (lT >= lt[0]) & (lT <= lt[-1])
        pcol = np.full(n_T, np.nan)
        ecol = np.full(n_T, np.nan)
        pcol[m] = 10.0 ** fp(lT[m])
        ecol[m] = fe(lT[m])
        partial.append((np.log10(r), m, pcol, ecol))
        # leave-one-out on interior points of this isochore. e error is
        # normalised by the isochore's e-SPAN, not |e| — e crosses zero
        # (binding-energy zero point), where a relative error is ill-defined.
        espan = Es.max() - Es.min()
        for k in range(1, len(Ts) - 1):
            sel = np.arange(len(Ts)) != k
            try:
                pk = 10.0 ** PchipInterpolator(lt[sel], np.log10(Ps[sel]))(lt[k])
                ek = PchipInterpolator(lt[sel], Es[sel])(lt[k])
            except ValueError:
                continue
            loo_err.append((abs(pk / Ps[k] - 1.0),
                            abs(ek - Es[k]) / max(espan, 1e-30)))

    # pass b: across isochores at each target T
    p = np.full((n_rho, n_T), np.nan)
    e = np.full((n_rho, n_T), np.nan)
    for j in range(n_T):
        lr = np.array([q[0] for q in partial if q[1][j]])
        if len(lr) < 4:
            continue
        o = np.argsort(lr)
        pv = np.array([q[2][j] for q in partial if q[1][j]])[o]
        ev = np.array([q[3][j] for q in partial if q[1][j]])[o]
        lr = lr[o]
        m = (lrho >= lr[0]) & (lrho <= lr[-1])
        p[m, j] = 10.0 ** PchipInterpolator(lr, np.log10(pv))(lrho[m])
        e[m, j] = PchipInterpolator(lr, ev)(lrho[m])

    hull = np.isfinite(p) & np.isfinite(e)
    # nearest-hull fill: along T within each row, then along rho per column
    for arr in (p, e):
        for i in range(n_rho):
            row = arr[i]
            good = np.where(np.isfinite(row))[0]
            if len(good):
                row[: good[0]] = row[good[0]]
                row[good[-1] + 1:] = row[good[-1]]
        for j in range(n_T):
            col = arr[:, j]
            good = np.where(np.isfinite(col))[0]
            if len(good):
                col[: good[0]] = col[good[0]]
                col[good[-1] + 1:] = col[good[-1]]
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(e)), "fill failed"

    loo = np.array(loo_err) if loo_err else np.zeros((0, 2))
    return lrho, lT, p, e, hull.astype(float), loo


def condition(lrho, lT, p, e, cv_floor_frac=1e-3):
    """Derivative blocks from PCHIP slopes on the conditioned surfaces.

    cv floored at cv_floor_frac * (3/2 kB / m_D) (ideal-gas fraction);
    dpdrho monotonised >= tiny positive. Counts recorded for provenance.
    Derivatives are wrt LINEAR rho/T (chain rule from the log axes).
    """
    rho = 10.0 ** lrho
    T = 10.0 ** lT
    nr, nt = p.shape
    dpdT = np.empty_like(p)
    cv = np.empty_like(p)
    dpdrho = np.empty_like(p)
    dedrho = np.empty_like(p)
    ln10 = np.log(10.0)
    for i in range(nr):
        dpdT[i] = PchipInterpolator(lT, p[i]).derivative()(lT) / (T * ln10)
        cv[i] = PchipInterpolator(lT, e[i]).derivative()(lT) / (T * ln10)
    for j in range(nt):
        dpdrho[:, j] = PchipInterpolator(lrho, p[:, j]).derivative()(lrho) / (rho * ln10)
        dedrho[:, j] = PchipInterpolator(lrho, e[:, j]).derivative()(lrho) / (rho * ln10)

    # Maxwell-loop check BEFORE monotonising (record only; PIMC tables are
    # expected loop-free — implement construction only if this trips)
    n_loops = int(np.sum(dpdrho < 0))

    cv_min = cv_floor_frac * 1.5 * KB / M_D
    n_cv = int(np.sum(cv < cv_min))
    cv = np.maximum(cv, cv_min)
    tiny = 1e-30
    n_mono = int(np.sum(dpdrho < tiny))
    dpdrho = np.maximum(dpdrho, tiny)
    # dpdT >= 0 is thermodynamically typical but not guaranteed; leave it.
    stats = dict(cv_floor=cv_min, cv_floored=n_cv, monotonised=n_mono,
                 maxwell=("none-needed" if n_loops == 0 else "LOOPS=%d" % n_loops))
    return dpdT, dpdrho, cv, dedrho, stats


def inverse_maps(lrho, lT, p, e, n_e=None, n_p=None):
    """Pre-inverted seed maps T(rho,e), T(rho,p) (Athena++-informed, D5/D8).

    Seed-only quality: each row's e(T)/p(T) is made monotone by cummax
    before inverting (documented; the runtime Newton polish against the
    forward surface supplies the accuracy and consistency).
    """
    T = 10.0 ** lT
    n_e = n_e or len(lT)
    n_p = n_p or len(lT)
    le = np.linspace(np.log10(e.min()), np.log10(e.max()), n_e)
    lp = np.linspace(np.log10(p.min()), np.log10(p.max()), n_p)
    T_of_e = np.empty((len(lrho), n_e))
    T_of_p = np.empty((len(lrho), n_p))
    for i in range(len(lrho)):
        em = np.maximum.accumulate(e[i])
        pm = np.maximum.accumulate(p[i])
        em += np.arange(len(em)) * 1e-12 * max(abs(em[-1]), 1.0)  # strictify
        pm *= 1.0 + np.arange(len(pm)) * 1e-12
        T_of_e[i] = np.interp(10.0 ** le, em, T, left=T[0], right=T[-1])
        T_of_p[i] = np.interp(10.0 ** lp, pm, T, left=T[0], right=T[-1])
    return le, lp, T_of_e, T_of_p


def shift_energy(e, hull):
    """Return (shifted e, shift) so min(e) on the hull is safely > 0."""
    emin = e[hull > 0.5].min()
    span = e[hull > 0.5].max() - emin
    shift = -emin + 1e-3 * span if emin <= 0 else 0.0
    return e + shift, shift


# =============================================================================
# QA
# =============================================================================


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
        # Parameterise the locus by T (its natural parameter): at each grid
        # temperature, solve for the density where the Rankine-Hugoniot
        # energy relation closes. Iterating over grid *densities* instead
        # undersamples badly — the whole in-table locus spans only ~0.1 dex
        # in rho (compression saturates near 4.7). Points whose Hugoniot T
        # lies below the table floor (the low-P foot) are absent by
        # construction.
        hr, hp = [], []
        for j in range(len(lT)):
            fe = PchipInterpolator(lrho, e[:, j])
            fp = PchipInterpolator(lrho, p[:, j])

            def resid(lr):
                return (fe(lr) - e0) - 0.5 * (fp(lr) + p0) * (
                    1.0 / rho0 - 1.0 / 10.0 ** lr)

            lo_lim = np.log10(rho0 * 1.02)
            m = lrho > lo_lim
            if m.sum() < 2:
                continue
            lrs = lrho[m]
            rvals = resid(lrs)
            for k in np.where(rvals[:-1] * rvals[1:] <= 0)[0]:
                lo, hi = lrs[k], lrs[k + 1]
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if resid(lo) * resid(mid) <= 0:
                        hi = mid
                    else:
                        lo = mid
                root = 0.5 * (lo + hi)
                hr.append(10.0 ** root / rho0)
                hp.append(float(fp(root)))
        o = np.argsort(hp)
        hr, hp = list(np.array(hr)[o]), list(np.array(hp)[o])
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


# =============================================================================
# commands
# =============================================================================


def provenance(material, source, comp, e_shift, stats):
    now = datetime.date.today().isoformat()
    return [
        ("material", material),
        ("source", source),
        ("generator", "eos_table_prep.py %s, run %s" % (git_sha(), now)),
        ("composition", comp),
        ("units", "cgs"),
        ("e_shift", "%.10e" % e_shift),
        ("conditioning", "cv_floor=%.4e cv_floored=%d monotonised=%d maxwell=%s"
         % (stats["cv_floor"], stats["cv_floored"], stats["monotonised"],
            stats["maxwell"])),
    ]


def cmd_synthetic(args):
    """Ideal-gas gamma-law table with analytic blocks (closed-form truth)."""
    g = args.gamma
    R = KB / M_D
    lrho = np.linspace(np.log10(args.rho_min), np.log10(args.rho_max), args.n_rho)
    lT = np.linspace(np.log10(args.T_min), np.log10(args.T_max), args.n_T)
    rho = 10.0 ** lrho[:, None]
    T = 10.0 ** lT[None, :]
    e = R * T / (g - 1.0) * np.ones_like(rho)
    p = rho * R * T
    blocks = {
        "p": p, "e": e,
        "dpdT": rho * R * np.ones_like(T),
        "dpdrho": R * T * np.ones_like(rho),
        "cv": R / (g - 1.0) * np.ones_like(p),
        "dedrho": np.zeros_like(p),
        "hull": np.ones_like(p),
    }
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks["T_of_e"], blocks["T_of_p"] = T_of_e, T_of_p
    stats = dict(cv_floor=0.0, cv_floored=0, monotonised=0, maxwell="none-needed")
    prov = provenance("ideal-D", "synthetic gamma-law gamma=%g (this tool)" % g,
                      "A=2.014 Z=1", 0.0, stats)
    write_eostab(args.out, prov, lrho, lT, blocks, {"le": le, "lp": lp})
    if args.qa:
        qa_plots(args.qa, lrho, lT, p, e, blocks["cv"], blocks["hull"],
                 tag="_synthetic")


def cmd_fpeos(args):
    raw = read_fpeos(args.src)
    pts = fpeos_to_deuterium(raw)
    print("source points: %d  (rho %.4g..%.4g g/cc D-equivalent, T %.4g..%.4g K)"
          % (len(pts["rho"]), pts["rho"].min(), pts["rho"].max(),
             pts["T"].min(), pts["T"].max()))
    lrho, lT, p, e, hull, loo = regrid(pts, args.n_rho, args.n_T)
    e, e_shift = shift_energy(e, hull)
    pts_plot = dict(pts, e=pts["e"])  # raw (unshifted) for overlays
    dpdT, dpdrho, cv, dedrho, stats = condition(lrho, lT, p, e)
    le, lp, T_of_e, T_of_p = inverse_maps(lrho, lT, p, e)
    blocks = dict(p=p, e=e, dpdT=dpdT, dpdrho=dpdrho, cv=cv, dedrho=dedrho,
                  T_of_e=T_of_e, T_of_p=T_of_p, hull=hull)
    prov = provenance(
        "D",
        "FPEOS %s (militzer.berkeley.edu, PRE 103 013203 (2021); H table "
        "isotope-scaled to D, see data/raw/README.md)" % os.path.basename(args.src),
        "A=2.014 Z=1", e_shift, stats)
    write_eostab(args.out, prov, lrho, lT, blocks, {"le": le, "lp": lp})
    print("conditioning: %s ; hull coverage %.1f%%"
          % (stats, 100.0 * hull.mean()))
    if args.qa:
        hug = dict(rho0=pts["hug_rho0"], e0=pts["hug_e0"] + e_shift,
                   p0=pts["hug_p0"])
        if args.hug_ref and os.path.exists(args.hug_ref):
            hug["ref"] = np.loadtxt(args.hug_ref)
            hug["ref_label"] = os.path.basename(args.hug_ref)
        qa_plots(args.qa, lrho, lT, p, e, cv, hull, pts=pts_plot, hug=hug,
                 loo=loo, tag="_D_fpeos")


def cmd_qa(args):
    prov, meta, blocks = read_eostab(args.table)
    qa_plots(args.qa, meta["lrho"], meta["lT"], blocks["p"], blocks["e"],
             blocks["cv"], blocks["hull"], tag="_" + os.path.basename(args.table))
    print("provenance:", {k: prov[k] for k in
                          ("material", "source", "conditioning") if k in prov})


def main():
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("synthetic", help="analytic ideal-gas .eostab")
    s.add_argument("--out", required=True)
    s.add_argument("--qa", default=None)
    s.add_argument("--gamma", type=float, default=1.4)
    s.add_argument("--n-rho", type=int, default=64)
    s.add_argument("--n-T", type=int, default=64)
    s.add_argument("--rho-min", type=float, default=1e-4)
    s.add_argument("--rho-max", type=float, default=10.0)
    s.add_argument("--T-min", type=float, default=1e3)
    s.add_argument("--T-max", type=float, default=1e7)
    s.set_defaults(func=cmd_synthetic)

    s = sub.add_parser("fpeos", help="condition an FPEOS element table -> D")
    s.add_argument("--src", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--qa", default=None)
    s.add_argument("--n-rho", type=int, default=96)
    s.add_argument("--n-T", type=int, default=96)
    s.add_argument("--hug-ref", default=None,
                   help="published Hugoniot points file (T rho P_GPa "
                        "compression) to overlay")
    s.set_defaults(func=cmd_fpeos)

    s = sub.add_parser("qa", help="QA plots for an existing .eostab")
    s.add_argument("--table", required=True)
    s.add_argument("--qa", required=True)
    s.set_defaults(func=cmd_qa)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
