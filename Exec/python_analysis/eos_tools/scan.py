"""G1 surface-monotonicity scan (doc/eos_table_conditioning_requirements.md).

The identified QA gap: the convexity gate checks cs^2 from the *floored
derivative blocks* and provably misses folds in the p value surface. This
module sweeps the surface itself — every isotherm in rho, every isochore
in T — and classifies each isotherm for the conditioning router:

    'monotone' : strictly increasing (no conditioning needed)
    'flat'     : no decreases, but exact/near tie-line flats (Maxwell-
                 constructed 311-style data) -> crossover replacement
    'loop'     : decreasing intervals (van der Waals loops, fill
                 collisions) -> Maxwell construction, then crossover

Run twice per table: on the raw source (routing) and on the emitted table
(acceptance gate: zero decreasing intervals).
"""

import numpy as np


def _rel_diff(p, axis):
    """(dp, scale) along axis with a sign-safe relative scale."""
    dp = np.diff(p, axis=axis)
    a = np.abs(p)
    sl_lo = [slice(None)] * p.ndim
    sl_hi = [slice(None)] * p.ndim
    sl_lo[axis] = slice(None, -1)
    sl_hi[axis] = slice(1, None)
    scale = np.maximum(np.maximum(a[tuple(sl_lo)], a[tuple(sl_hi)]), 1e-300)
    return dp, scale


def classify_isotherm(p, flat_tol=1e-13):
    """Classify one isotherm p(rho) -> 'monotone' | 'flat' | 'loop'."""
    dp = np.diff(p)
    scale = np.maximum(np.maximum(np.abs(p[:-1]), np.abs(p[1:])), 1e-300)
    if np.any(dp < -flat_tol * scale):
        return "loop"
    if np.any(np.abs(dp) <= flat_tol * scale) or np.any(p <= 0.0):
        return "flat"  # nonpositive cells are treated as band cells too
    return "monotone"


def g1_scan(p, flat_tol=1e-13):
    """Scan a (Nr, Nt) pressure surface. Returns a stats dict.

    Gate (emitted table): ok_rho and ok_T — zero decreasing intervals in
    either direction. worst_* are the largest relative decreases found.
    """
    dp_r, sc_r = _rel_diff(p, axis=0)
    dp_t, sc_t = _rel_diff(p, axis=1)
    dec_r = dp_r < -flat_tol * sc_r
    dec_t = dp_t < -flat_tol * sc_t
    per = [classify_isotherm(p[:, j], flat_tol) for j in range(p.shape[1])]
    return {
        "per_isotherm": per,
        "n_loop": per.count("loop"),
        "n_flat": per.count("flat"),
        "n_monotone": per.count("monotone"),
        "n_dec_rho": int(dec_r.sum()),
        "n_dec_T": int(dec_t.sum()),
        "worst_drop_rho": float((-dp_r / sc_r).max(initial=0.0)),
        "worst_drop_T": float((-dp_t / sc_t).max(initial=0.0)),
        "n_nonpos": int(np.sum(p <= 0.0)),
        "ok_rho": not dec_r.any(),
        "ok_T": not dec_t.any(),
    }


def report(stats, label=""):
    print("G1[%s] isotherms: %d monotone / %d flat / %d loop; "
          "dec intervals rho=%d (worst %.2e) T=%d (worst %.2e); nonpos=%d"
          % (label, stats["n_monotone"], stats["n_flat"], stats["n_loop"],
             stats["n_dec_rho"], stats["worst_drop_rho"],
             stats["n_dec_T"], stats["worst_drop_T"], stats["n_nonpos"]))
