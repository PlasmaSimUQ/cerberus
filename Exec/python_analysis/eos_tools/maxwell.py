"""Maxwell construction + steep monotone crossover (H1/H2, UC1-d).

Replaces every two-phase feature of a p(rho)|_T isotherm — van der Waals
loops, exact Maxwell tie-line flats, nonpositive (tension) cells — with a
strictly monotone-in-rho bridge, so that rho(p,T) and T(rho,p) are
single-rooted everywhere (requirement H1) and no pressure is ever pinned
flat (H2, G4: tension_clip = 0 by construction).

Route per isotherm (doc/eos_sesame_plan.md §3, stage 3-4):
  1. 'loop'  -> equal-area Maxwell construction locates the tie line
               (p*, band); where the construction is numerically
               degenerate (vapor pressure underflow at deep-cold T, no
               admissible p* > 0) the band is bridged directly.
  2. 'flat'  -> the shipped tie line (311-style) is the band.
  3. bridge  -> band cells are replaced by a strictly increasing
               geometric ramp between the branch anchor values; e across
               any replaced band follows the two-phase lever rule
               (linear in specific volume v = 1/rho between the band-edge
               values), preserving the latent-heat energy content.

The decided cavitated-response stiffness (plan §3.0, D-e) is applied to
the DERIVATIVE BLOCKS by the caller, never to this value surface.
"""

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy 2.x / 1.x


def equal_area_pstar(rho, p):
    """Equal-area Maxwell pressure for one folded isotherm.

    rho ascending; p carries a van der Waals loop (some d(p)/d(rho) < 0).
    Works in v = 1/rho. Returns (pstar, i_lo, i_hi) with [i_lo, i_hi] the
    inclusive rho-index band strictly inside the tie line, or None when
    the construction is degenerate (no admissible pstar > 0, fewer than
    two crossings, or the area residual cannot be bracketed).
    """
    v = 1.0 / rho[::-1]          # ascending v (gas to the right)
    pv = p[::-1].astype(float)
    n = len(v)
    dpdv = np.diff(pv)
    up = np.where(dpdv > 0.0)[0]  # unstable/loop segments: p rising with v
    if len(up) == 0:
        return None
    # tie-line pressure bracket: the loop's local extrema
    p_hi = pv[up[0]:up[-1] + 2].max()
    p_lo = pv[up[0]:up[-1] + 2].min()
    p_lo = max(p_lo, p_hi * 1e-12, 0.0)
    if not (p_hi > 0.0) or p_lo >= p_hi:
        return None

    def crossings(P):
        s = pv - P
        idx = np.where(s[:-1] * s[1:] < 0.0)[0]
        xs = [v[i] + (v[i + 1] - v[i]) * s[i] / (s[i] - s[i + 1]) for i in idx]
        for i in np.where(s == 0.0)[0]:
            xs.append(v[i])
        return sorted(xs)

    def area(P):
        xs = crossings(P)
        if len(xs) < 2:
            return None
        va, vb = xs[0], xs[-1]
        m = (v > va) & (v < vb)
        vv = np.concatenate(([va], v[m], [vb]))
        pp = np.concatenate(([P], pv[m], [P]))
        return _trapz(pp - P, vv)

    lo, hi = p_lo * (1.0 + 1e-9), p_hi * (1.0 - 1e-9)
    a_lo, a_hi = area(lo), area(hi)
    if a_lo is None or a_hi is None or a_lo * a_hi > 0.0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        a = area(mid)
        if a is None:
            return None
        if a_lo * a > 0.0:
            lo, a_lo = mid, a
        else:
            hi = mid
        if (hi - lo) <= 1e-14 * hi:
            break
    pstar = 0.5 * (lo + hi)
    xs = crossings(pstar)
    if len(xs) < 2:
        return None
    va, vb = xs[0], xs[-1]
    # v strictly inside (va, vb) -> rho strictly inside (1/vb, 1/va)
    inside = np.where((rho > 1.0 / vb) & (rho < 1.0 / va))[0]
    if len(inside) == 0:
        return None
    return pstar, int(inside[0]), int(inside[-1])


def _geom_bridge(x, xa, ya, xb, yb):
    """Strictly increasing bridge y(x) between (xa, ya) and (xb, yb),
    geometric (log-linear) when both anchors are positive."""
    t = (x - xa) / (xb - xa)
    if ya > 0.0 and yb > 0.0 and yb > ya:
        return ya * (yb / ya) ** t
    return ya + (yb - ya) * t


KNEE_GAP = 1.5  # min anchor/foot ratio for the hug shape to engage


def crossover_isotherm(rho, p, e, flat_tol=1e-13, rel_eps=1e-12,
                       p_foot=None, allow_maxwell=True):
    """Condition one isotherm: strictly-monotone p(rho), lever-rule e.

    Returns (p2, e2, mask, info): mask flags every cell whose value was
    replaced (-> hull 0 in the emitted table); info records the route
    ('none' | 'flat' | 'maxwell' | 'ramp') and counts.

    p_foot (cgs, quiet-foot plan A): when set, every rebuilt bridge whose
    left edge is below p_foot and whose anchor exceeds KNEE_GAP*p_foot
    hugs p_foot at the top of the replaced span (the knee), confining the
    steep rise to the single cell before the anchor — instead of the
    legacy log-linear ramp across the whole span, which parks
    anchor-scale (GPa) values at the solid's own density. Legacy shape is
    byte-identical when p_foot is None.

    allow_maxwell=False routes folded isotherms to the ramp rebuild
    without an equal-area construction — required for the post-blend
    crossover in coldext, where negatives are construction tension and a
    tie line is meaningless (quiet-foot plan B / Al 96.5 GPa defect).
    """
    lrho = np.log10(rho)
    p2 = p.astype(float).copy()
    e2 = e.astype(float).copy()
    n = len(p2)
    info = {"route": "none", "n_replaced": 0, "pstar": None}

    dp = np.diff(p2)
    scale = np.maximum(np.maximum(np.abs(p2[:-1]), np.abs(p2[1:])), 1e-300)
    has_dec = np.any(dp < -flat_tol * scale)
    has_flat = np.any(np.abs(dp) <= flat_tol * scale)
    has_nonpos = np.any(p2 <= 0.0)
    if not (has_dec or has_flat or has_nonpos):
        return p2, e2, np.zeros(n, bool), info

    # 1. Maxwell placement on folded isotherms (before any repair, while
    #    the loop is exactly representable)
    if has_dec:
        mx = (equal_area_pstar(rho, p2)
              if allow_maxwell and not has_nonpos else None)
        if mx is not None:
            pstar, i0, i1 = mx
            p2[i0:i1 + 1] = pstar
            info.update(route="maxwell", pstar=pstar)
        else:
            info["route"] = "ramp"
    else:
        info["route"] = "flat"

    # 2. nonpositive (tension) cells: graded tiny-positive placeholders so
    #    the monotone rebuild below has admissible values to bridge; every
    #    such cell is band (mask) by definition
    pos = p2 > 0.0
    if not pos.any():
        raise ValueError("isotherm with no positive pressure anywhere — "
                         "trim the T range (--T-min) below this isotherm")
    if has_nonpos:
        p_ref = p2[pos].min()
        bad = np.where(~pos)[0]
        # ascending tiny values far below the smallest real pressure
        p2[bad] = p_ref * 1e-10 * (1.0 + 1e-6 * np.arange(len(bad)))

    # 3. strictly-monotone rebuild: upper envelope, then geometric bridges
    #    across every run the envelope left flat
    pm = np.maximum.accumulate(p2)
    i = 0
    while i < n - 1:
        if pm[i + 1] > pm[i] * (1.0 + flat_tol):
            i += 1
            continue
        j = i + 1
        while j < n - 1 and pm[j + 1] <= pm[i] * (1.0 + flat_tol):
            j += 1
        if j < n - 1:
            if (p_foot is not None and pm[i] < p_foot
                    and pm[j + 1] > KNEE_GAP * p_foot):
                # floor-hugging shape: knee at the top of the span holds
                # p_foot; the rise to the anchor is the final cell only
                pm[i + 1:j + 1] = _geom_bridge(lrho[i + 1:j + 1], lrho[i],
                                               pm[i], lrho[j], p_foot)
            else:
                pm[i + 1:j + 1] = _geom_bridge(lrho[i + 1:j + 1], lrho[i],
                                               pm[i], lrho[j + 1], pm[j + 1])
        else:  # flat run reaches the top of the grid: epsilon up-ramp
            pm[i + 1:] = pm[i] * (1.0 + rel_eps) ** np.arange(1, n - i)
        i = j
    mask = np.abs(pm - p) > 10.0 * flat_tol * np.maximum(np.abs(p), 1e-300)
    mask |= p <= 0.0
    p2 = pm

    # 4. lever-rule e across every replaced band (linear in v between the
    #    band-edge anchors — the two-phase mixture energy, latent heat kept)
    k = 0
    while k < n:
        if not mask[k]:
            k += 1
            continue
        m0 = k
        while k < n and mask[k]:
            k += 1
        m1 = k - 1
        a, b = m0 - 1, m1 + 1
        if a >= 0 and b < n and m1 > m0:
            v = 1.0 / rho
            t = (v[m0:m1 + 1] - v[a]) / (v[b] - v[a])
            e2[m0:m1 + 1] = e[a] + t * (e[b] - e[a])
    info["n_replaced"] = int(mask.sum())
    return p2, e2, mask, info


P_FOOT_TILT = 1e-3  # strict-dpdT tilt on the hugged foot (quiet-foot D1)


def condition_surface(rho, T, p, e, flat_tol=1e-13, p_foot=None,
                      allow_maxwell=True):
    """Apply crossover_isotherm to every column of p, e (Nr, Nt).

    Returns (p2, e2, mask, stats); stats feed the conditioning provenance
    line (routes, replaced-cell count, tie-line pressures range).

    p_foot (cgs): quiet-foot target, tilted upward by P_FOOT_TILT across
    the T axis so the hugged foot stays strictly increasing in T (no
    exact flats for the (rho, p) -> T inversion to trip on).
    """
    nr, nt = p.shape
    p2 = np.empty_like(p, dtype=float)
    e2 = np.empty_like(e, dtype=float)
    mask = np.zeros((nr, nt), bool)
    routes = {"none": 0, "flat": 0, "maxwell": 0, "ramp": 0}
    pstars = []
    for j in range(nt):
        pf = (p_foot * (1.0 + P_FOOT_TILT * j / max(nt - 1, 1))
              if p_foot is not None else None)
        pj, ej, mj, info = crossover_isotherm(rho, p[:, j], e[:, j], flat_tol,
                                              p_foot=pf,
                                              allow_maxwell=allow_maxwell)
        p2[:, j], e2[:, j], mask[:, j] = pj, ej, mj
        routes[info["route"]] += 1
        if info["pstar"] is not None:
            pstars.append(info["pstar"])
    stats = dict(routes=routes, n_band_cells=int(mask.sum()),
                 n_maxwell=routes["maxwell"], n_ramp=routes["ramp"],
                 n_flat=routes["flat"],
                 pstar_min=(min(pstars) if pstars else None),
                 pstar_max=(max(pstars) if pstars else None))
    return p2, e2, mask, stats
