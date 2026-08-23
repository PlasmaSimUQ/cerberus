"""Sub-floor temperature extension — rows below a table's native T floor.

Purpose: a mixture is bracketed in T by the HIGHEST native floor among its
retained components; this stage lowers one member's floor (first use: the
eref295 set member with a 100 K native floor, brought to the set's ~17.8 K).
Plan: doc/eos_air_lowT_extension_plan.md. Decisions as built (plan labels):
D1-a (rows prepended on the NATIVE log spacing, native nodes untouched),
D2-a (uniform ideal-thermal branch), D3-a (new rows are hull 0 —
constructed, not data), D6-a (first-derivative seam kink accepted).

Construction, per density column i with native floor row (T0, p0, e0):

    e(rho, T) = e0 - cv_x * (T0 - T)                  cv_x = 5/2 k/m
    p(rho, T) = max( p0 - rho*R_s*(T0 - T),           anchored thermal slope
                     p0 * T/T0 )                      scaled-row floor
                                                      R_s = k/m, m molecular

Why this pair (all provable, none tuned):
  * vapor columns (p0 = rho R_s T0): both branches equal rho R_s T — the
    ideal gas exactly;
  * condensed columns (p0 > rho R_s T0): anchored branch wins, slope
    rho R_s, p stays dominated by the cold curve in p0;
  * two-phase columns of the anchor row (p0 = p_sat < rho R_s T0): the
    scaled branch wins — the ideal-gas law at the dome-edge vapor density,
    an upper bound on p_sat(T), positive and increasing in T;
  * dp/dT = rho R_s on branch 1 and p0/T0 on branch 2 — equal exactly at
    the switch (p0 = rho R_s T0), so dpdT is continuous and > 0 everywhere;
  * both branches are <= p0 for T < T0, so the cummax-in-T enforcement can
    never lift a native row (native byte identity is structural, not gated);
  * dp/drho >= 0 on both branches wherever p0 is non-decreasing in rho (it
    is, post-monotonise); strictness on flat (dome) spans is restored by
    the existing monotonise_rho epsilon ramps, as for the anchor row itself.

Limitations are the plan's §6: below the anchor row everything is
construction — metastable ideal vapor (no condensation/latent heat), gas
cv and Gamma at condensed density, no solid phases.
"""
import numpy as np

from .constants import AMU_G, KB


def extend_T_floor(lrho, lT, p, e, hull, band, lT_floor, m_amu):
    """Prepend rows below lT[0] down to (at or just below) lT_floor.

    lrho, lT: log10 axes (lT ascending, uniform). p, e, hull, band:
    (Nr, Nt) arrays in cgs. m_amu: molecular mass (amu) for R_s and cv_x.

    Returns (lT2, p2, e2, hull2, band2, st). New rows: hull 0, band False.
    st: rows (k), lT_floor (achieved), h (native spacing), R_s, cv_x,
    n_vapor / n_anchor / n_scaled (new cells per regime), p_min (new rows).
    A no-op (k = 0) when the native floor is already at/below lT_floor.
    """
    lrho = np.asarray(lrho, dtype=float)
    lT = np.asarray(lT, dtype=float)
    if len(lT) < 2:
        raise ValueError("extend_T_floor needs at least two T nodes")
    h = float(lT[1] - lT[0])
    if h <= 0:
        raise ValueError("lT must be ascending")
    st = dict(rows=0, lT_floor=float(lT[0]), h=h,
              R_s=KB / (m_amu * AMU_G), cv_x=2.5 * KB / (m_amu * AMU_G),
              n_vapor=0, n_anchor=0, n_scaled=0, p_min=float("nan"))
    k = int(np.ceil((lT[0] - lT_floor) / h - 1e-9))
    if k <= 0:
        return lT, p, e, hull, band, st
    lT_new = lT[0] - h * np.arange(k, 0, -1)          # ascending, below lT[0]
    T0 = 10.0 ** lT[0]
    T = 10.0 ** lT_new
    rho = 10.0 ** lrho
    R_s, cv_x = st["R_s"], st["cv_x"]

    p0 = p[:, 0][:, None]
    e0 = e[:, 0][:, None]
    dT = (T0 - T)[None, :]
    p_anchor = p0 - rho[:, None] * R_s * dT
    p_scaled = p0 * (T / T0)[None, :]
    p_new = np.maximum(p_anchor, p_scaled)
    e_new = e0 - cv_x * dT

    p2 = np.concatenate([p_new, p], axis=1)
    e2 = np.concatenate([e_new, e], axis=1)
    hull2 = np.concatenate([np.zeros_like(p_new), hull], axis=1)
    band2 = np.concatenate([np.zeros(p_new.shape, dtype=bool), band], axis=1)
    lT2 = np.concatenate([lT_new, lT])
    # regime tally per cell: vapor (both branches equal — the ideal gas),
    # anchored (condensed, p0 > rho R_s T0), scaled (sub-ideal anchor row)
    ideal0 = rho[:, None] * R_s * T0
    vapor = np.abs(p0 - ideal0) <= 1e-9 * ideal0
    vapor = np.broadcast_to(vapor, p_new.shape)
    st.update(rows=k, lT_floor=float(lT_new[0]),
              n_vapor=int(np.sum(vapor)),
              n_anchor=int(np.sum(~vapor & (p_anchor >= p_scaled))),
              n_scaled=int(np.sum(~vapor & (p_anchor < p_scaled))),
              p_min=float(p_new.min()))
    return lT2, p2, e2, hull2, band2, st
