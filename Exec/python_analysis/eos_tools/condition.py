"""Conditioning of regridded (p, e) surfaces: derivative blocks, inverse
maps, energy shift. All curation happens here, offline (plan D4)."""

import numpy as np
from scipy.interpolate import PchipInterpolator

from .constants import KB, M_D


def monotonise_T(F, rel_eps=1e-12):
    """Enforce strictly nondecreasing F(T) along every rho-row.

    The v1 inversion contract (single-branch invert_T_from_e/p) needs
    monotone value surfaces; the conditioning already floors cv > 0, so
    this makes the values consistent with the derivative blocks. Measured
    magnitudes are small (worst p dip 1.2e-2 relative, e 3.6e-4, in the
    REOS.3 near-melt rows). Returns (F', n_changed, max_rel_change).
    (Moved here from splice.py so surface conditioning needs no splice
    imports; splice.py re-exports it.)
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


def monotonise_rho(F, rel_eps=1e-12):
    """Enforce strictly nondecreasing F(rho) along every isotherm (the H1
    surface contract, in-ρ direction). Same construction as monotonise_T."""
    Fm, n, mrel = monotonise_T(np.ascontiguousarray(F.T), rel_eps)
    return np.ascontiguousarray(Fm.T), n, mrel


def condition(lrho, lT, p, e, cv_floor_frac=1e-3, m_ref=M_D,
              cs_floor_cms=None):
    """Derivative blocks from PCHIP slopes on the conditioned surfaces.

    cv floored at cv_floor_frac * (3/2 kB / m_ref) (ideal-gas fraction;
    m_ref defaults to the deuterium mass for the existing D pipelines —
    pass the material's mean particle mass for other materials);
    dpdrho monotonised >= tiny positive. Counts recorded for provenance.
    Derivatives are wrt LINEAR rho/T (chain rule from the log axes).

    cs_floor_cms (W-B, ti_splice_plan_v2 §0.3): OPT-IN sound-speed safety
    net in cm/s (never code units — those are refs-dependent). Implemented
    as dpdrho >= cs_floor^2; since cs^2 = dpdrho + T*dpdT^2/(rho^2 cv) and
    the added term is >= 0 once cv is floored positive, this guarantees
    cs >= cs_floor everywhere. Off by default: a default-on floor would
    silently break the trackP byte-identity baseline (T1).
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

    cv_min = cv_floor_frac * 1.5 * KB / m_ref
    n_cv = int(np.sum(cv < cv_min))
    cv = np.maximum(cv, cv_min)
    tiny = 1e-30
    n_mono = int(np.sum(dpdrho < tiny))
    dpdrho = np.maximum(dpdrho, tiny)
    n_cs = 0
    if cs_floor_cms:
        c2 = float(cs_floor_cms) ** 2
        n_cs = int(np.sum(dpdrho < c2))
        dpdrho = np.maximum(dpdrho, c2)
    # dpdT >= 0 is thermodynamically typical but not guaranteed; leave it.
    stats = dict(cv_floor=cv_min, cv_floored=n_cv, monotonised=n_mono,
                 cs_floor=(float(cs_floor_cms) if cs_floor_cms else 0.0),
                 cs_floored=n_cs,
                 maxwell=("none-needed" if n_loops == 0 else "LOOPS=%d" % n_loops))
    return dpdT, dpdrho, cv, dedrho, stats


def thermal_floor(dpdrho, lT, m_ref):
    """Ideal-gas isothermal stiffness floor (W-B): no single-phase fluid is
    isothermally softer than ideal gas, dpdrho >= kB*T/m.

    Returns (dpdrho_floored, kTm, sub_mask). SCOPE POLICY IS THE CALLER'S:
    the SESAME path applies the test to all cells and demotes failures to
    band/hull-0 (correct there — failures are tie-line residue); a splice
    path with in-hull data softer than the ATOMIC bound (e.g. molecular D2)
    must restrict the test to band/hull-0 cells (ti_splice_plan_v2 B7).
    """
    lT = np.asarray(lT, float)
    kTm = KB * (10.0 ** lT)[None, :] / m_ref * np.ones((dpdrho.shape[0], 1))
    sub = dpdrho < kTm
    return np.maximum(dpdrho, kTm), kTm, sub


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
