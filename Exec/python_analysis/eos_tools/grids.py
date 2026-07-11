"""Hull-aware regrid of scattered (rho, T) source points onto log-uniform grids."""

import numpy as np
from scipy.interpolate import PchipInterpolator


def isochores(pts):
    """Group scattered points into per-density isochores, T-sorted."""
    out = []
    for r in np.unique(pts["rho"]):
        m = pts["rho"] == r
        o = np.argsort(pts["T"][m])
        out.append((r, pts["T"][m][o], pts["p"][m][o], pts["e"][m][o]))
    return out


def regrid(pts, n_rho, n_T):
    """Hull-aware two-pass PCHIP regrid onto a log-uniform (rho, T) grid
    spanning the source's own range. See regrid_onto for the mechanics."""
    iso = isochores(pts)
    lrho = np.linspace(np.log10(min(r for r, *_ in iso)),
                       np.log10(max(r for r, *_ in iso)), n_rho)
    lT = np.linspace(np.log10(pts["T"].min()), np.log10(pts["T"].max()), n_T)
    p, e, hull, loo = regrid_onto(pts, lrho, lT)
    return lrho, lT, p, e, hull, loo


def regrid_onto(pts, lrho, lT):
    """Hull-aware two-pass PCHIP regrid onto GIVEN log10 axes.

    Pass a: along each source isochore in log T (within its own T range).
    Pass b: along log rho at each target T (within the covered rho range).
    Cells never covered are filled with the nearest hull value and marked
    0 in the hull mask. Refuses to extrapolate anywhere.
    Returns (p, e, hull, loo) with p,e shaped (len(lrho), len(lT));
    loo = leave-one-out regrid-error stats (review amendment: measure,
    don't assume, on a coarse source).
    """
    iso = isochores(pts)
    n_rho, n_T = len(lrho), len(lT)

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
    return p, e, hull.astype(float), loo
