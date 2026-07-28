"""G1 scan classification + Maxwell/crossover invariants.

The loop fixture is a van der Waals isotherm below critical: the
equal-area pressure is checked against an independent construction, and
the crossover output must satisfy exactly the properties the runtime
needs — strict monotonicity in rho, positivity, endpoint preservation,
lever-rule energy across the band.
"""

import numpy as np
import pytest

from eos_tools.maxwell import condition_surface, crossover_isotherm, equal_area_pstar
from eos_tools.scan import classify_isotherm, g1_scan

# van der Waals in reduced units: p = 8T/(3v-1) - 3/v^2, T=0.9 (subcritical)
TR = 0.9


def vdw_p(v):
    return 8.0 * TR / (3.0 * v - 1.0) - 3.0 / v ** 2


@pytest.fixture(scope="module")
def loop_isotherm():
    v = np.geomspace(0.55, 30.0, 400)
    rho = 1.0 / v[::-1]
    p = vdw_p(1.0 / rho)
    e = 3.0 * TR - 3.0 * rho  # vdW-like energy: e decreases with v via a=3
    return rho, p, e


def test_classify():
    assert classify_isotherm(np.array([1.0, 2.0, 3.0])) == "monotone"
    assert classify_isotherm(np.array([1.0, 2.0, 2.0, 3.0])) == "flat"
    assert classify_isotherm(np.array([1.0, 2.0, 1.5, 3.0])) == "loop"
    assert classify_isotherm(np.array([-1.0, 1.0, 2.0])) == "flat"  # nonpos


def test_equal_area_matches_reference(loop_isotherm):
    rho, p, _ = loop_isotherm
    res = equal_area_pstar(rho, p)
    assert res is not None
    pstar, i0, i1 = res
    # independent reference: dense scan of the area residual
    v = np.geomspace(0.55, 30.0, 20000)
    pv = vdw_p(v)

    def area(P):
        s = pv - P
        idx = np.where(s[:-1] * s[1:] < 0)[0]
        if len(idx) < 3:          # no tie line at this P: not a candidate
            return np.inf
        xs = [v[i] + (v[i + 1] - v[i]) * s[i] / (s[i] - s[i + 1]) for i in idx]
        va, vb = xs[0], xs[-1]
        m = (v > va) & (v < vb)
        return np.trapz(np.concatenate(([0.0], pv[m] - P, [0.0])),
                        np.concatenate(([va], v[m], [vb])))

    Ps = np.linspace(0.5 * pstar, 1.5 * pstar, 2001)
    ref = Ps[np.argmin(np.abs([area(P) for P in Ps]))]
    assert abs(pstar - ref) / ref < 2e-3
    assert 0 < i0 <= i1 < len(rho) - 1


def test_crossover_properties(loop_isotherm):
    rho, p, e = loop_isotherm
    p2, e2, mask, info = crossover_isotherm(rho, p, e)
    assert info["route"] == "maxwell"
    assert np.all(np.diff(p2) > 0)                    # strictly monotone
    assert np.all(p2 > 0)                             # H4-compatible
    assert p2[0] == p[0] and p2[-1] == p[-1]          # branches preserved
    assert e2[0] == e[0] and e2[-1] == e[-1]
    assert mask.any() and not mask[0] and not mask[-1]
    # lever rule: e linear in v across each replaced band
    band = np.where(mask)[0]
    a, b = band[0] - 1, band[-1] + 1
    v = 1.0 / rho
    t = (v[band] - v[a]) / (v[b] - v[a])
    assert np.allclose(e2[band], e[a] + t * (e[b] - e[a]), rtol=0, atol=1e-12)


def test_flat_and_tension_routes():
    rho = np.geomspace(1e-3, 10.0, 50)
    # shipped Maxwell flat (311-style): rising, exactly flat, rising
    p = np.where(rho < 0.1, rho, np.where(rho < 1.0, 0.1, 0.1 * rho))
    e = np.linspace(1.0, 2.0, 50)
    p2, e2, mask, info = crossover_isotherm(rho, p, e)
    assert info["route"] == "flat"
    assert np.all(np.diff(p2) > 0) and mask.any()
    # tension isotherm (301-style cold): negative pressures in the dome
    pt = np.where(rho < 1.0, -0.5 + 0.1 * rho, 2.0 * (rho - 1.0) + 1e-4)
    p2, e2, mask, info = crossover_isotherm(rho, pt, e)
    assert np.all(np.diff(p2) > 0) and np.all(p2 > 0)
    assert mask[0]                                    # nonpos cells banded
    with pytest.raises(ValueError, match="no positive pressure"):
        crossover_isotherm(rho, -np.ones(50), e)


def test_condition_surface_and_g1(loop_isotherm):
    rho, p, e = loop_isotherm
    T = np.array([0.9, 0.95, 1.2])
    # column 0: loop; column 1: same loop; column 2: monotone (supercritical)
    v = 1.0 / rho
    psup = 8.0 * 1.2 / (3.0 * v - 1.0) - 3.0 / v ** 2
    P = np.column_stack([p, p, np.sort(psup)])
    E = np.column_stack([e, e, e])
    raw = g1_scan(P)
    assert raw["n_loop"] == 2 and not raw["ok_rho"]
    p2, e2, mask, stats = condition_surface(rho, T, P, E)
    out = g1_scan(p2)
    assert out["ok_rho"] and out["n_dec_rho"] == 0    # the G1 gate
    assert stats["n_maxwell"] == 2 and stats["routes"]["none"] == 1
    assert not mask[:, 2].any()                       # monotone col untouched
