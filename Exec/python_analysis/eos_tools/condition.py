"""Conditioning of regridded (p, e) surfaces: derivative blocks, inverse
maps, energy shift. All curation happens here, offline (plan D4)."""

import numpy as np
from scipy.interpolate import PchipInterpolator

from .constants import KB, M_D


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
