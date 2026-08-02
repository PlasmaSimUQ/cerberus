"""cold_extend zone mechanics on a synthetic backbone (plan v2 T3/T4).

A toy 'SESAME' surface with a ragged hull edge and a toy solid model with
a known constant energy offset exercise: zone-1 replacement, the tanh
blend's seam-free decay, the gap-column bridge, the melt-cap exclusion of
vapour columns, and the H5 offset-constancy gate (recovery and failure).
"""

import numpy as np
import pytest

from eos_tools.coldext import N_HB, W_EPS, cold_extend, melt_curve_on

NR, NT = 24, 60
LRHO = np.linspace(0.0, 1.0, NR)          # rho 1..10 g/cc
LT = np.linspace(1.0, 3.0, NT)            # T 10..1000 K
RHO = 10.0 ** LRHO
T = 10.0 ** LT
E_OFF = 3.0e9                              # erg/g, the constant to recover


class ToyModel:
    m_atom = 1.0e-23

    def p(self, R, TT):
        return 2.0e8 * R * TT + 5.0e9

    def e(self, R, TT):
        return 4.0e7 * TT + 1.0e8 * R - E_OFF


def toy_sesame():
    """p, e smooth; hull edge at 40 K for rho >= 3 (overlap columns),
    at 600 K for rho < 3 (gap columns)."""
    R = RHO[:, None] * np.ones((1, NT))
    TT = T[None, :] * np.ones((NR, 1))
    p = 2.1e8 * R * TT + 5.5e9
    e = 4.0e7 * TT + 1.0e8 * R          # = model e + E_OFF exactly
    hull = np.zeros((NR, NT))
    for i in range(NR):
        edge = 40.0 if RHO[i] >= 3.0 else 600.0
        hull[i, T >= edge] = 1.0
    return p, e, hull


def melt():
    """f*Tm above the grid floor everywhere except the two lowest-density
    columns (vapour exclusion). Below rho = 3 the cap (270 K) sits under
    the 600 K hull edge -> gap columns with a bridge; above rho = 3 the
    cap (810 K) clears the 40 K edge + blend -> overlap columns."""
    Tm = np.where(RHO >= 3.0, 900.0, 300.0)
    Tm[:2] = 5.0                        # 0.9*5 K < 10 K floor -> untouched
    return Tm


def run(align_tol=5.0, model=None):
    p, e, hull = toy_sesame()
    band = np.zeros((NR, NT), bool)
    return cold_extend(LRHO, LT, p, e, hull, band, model or ToyModel(),
                       melt(), f_melt=0.9, align_tol=align_tol)


def test_alignment_recovers_constant_offset():
    _, _, _, _, st = run()
    assert st["align_n"] > 0
    assert st["align_cE"] == pytest.approx(E_OFF, rel=1e-12)
    assert st["align_std_kT"] < 1e-9


def test_zone1_is_aligned_model_and_hull0():
    p, e, hull, band, st = run()
    m = ToyModel()
    i = NR - 1                          # deep overlap column (rho = 10)
    jh = np.flatnonzero(toy_sesame()[2][i] > 0.5)[0]
    for j in range(jh):
        assert p[i, j] == pytest.approx(m.p(RHO[i], T[j]), rel=1e-12)
        assert e[i, j] == pytest.approx(m.e(RHO[i], T[j]) + E_OFF, rel=1e-12)
        assert hull[i, j] == 0.0
    assert st["zone1"] > 0 and st["blend"] > 0


def test_blend_decays_without_seam():
    p, _, hull, _, _ = run()
    p0 = toy_sesame()[0]
    i = NR - 1
    jh = np.flatnonzero(toy_sesame()[2][i] > 0.5)[0]
    # in-hull blend cells keep hull = 1
    assert np.all(hull[i, jh:] == 1.0)
    # at the last modified cell the weight is < W_EPS-ish: no step there
    changed = np.flatnonzero(p[i] != p0[i])
    jlast = changed[-1]
    rel_jump = abs(p[i, jlast] - p0[i, jlast]) / abs(p0[i, jlast])
    assert rel_jump < 10 * W_EPS
    # and the blend spans at least the R3 half-band above the edge
    assert jlast >= jh + N_HB


def test_gap_column_bridge_linear_and_hull0():
    p, e, hull, _, st = run()
    i = 5                               # rho ~ 1.6 < 3 -> gap column
    hull0 = toy_sesame()[2]
    jh = np.flatnonzero(hull0[i] > 0.5)[0]
    jc = np.searchsorted(T, 0.9 * 300.0, side="right") - 1
    assert jc < jh
    assert np.all(hull[i, :jh] == 0.0)
    # bridge is linear in log T between the two anchors
    tt = (LT[jc + 1:jh] - LT[jc]) / (LT[jh] - LT[jc])
    m = ToyModel()
    pa = m.p(RHO[i], T[jc])
    want = pa + tt * (p[i, jh] - pa)
    assert np.allclose(p[i, jc + 1:jh], want, rtol=1e-12)
    assert st["bridge"] > 0


def test_melt_cap_excludes_vapour_columns():
    p, e, hull, _, _ = run()
    p0, e0, hull0 = toy_sesame()
    for i in (0, 1):                    # f*Tm below the grid floor
        assert np.array_equal(p[i], p0[i])
        assert np.array_equal(e[i], e0[i])
        assert np.array_equal(hull[i], hull0[i])


def test_offset_constancy_gate_fails_on_T_dependence():
    class Drifting(ToyModel):
        def e(self, R, TT):
            return super().e(R, TT) - 5.0e7 * TT   # offset drifts with T

    with pytest.raises(RuntimeError, match="offset-constancy"):
        run(align_tol=0.5, model=Drifting())


def test_melt_curve_interpolation_clamps():
    ml = {"rho": np.array([0.0, 0.1, 1.0, 10.0]),
          "Tm": np.array([0.0, 10.0, 100.0, 1000.0])}
    Tm = melt_curve_on(np.array([-3.0, -1.0, 0.0, 1.0, 2.0]), ml)
    assert Tm[0] == pytest.approx(10.0)     # clamped below
    assert Tm[2] == pytest.approx(100.0)
    assert Tm[4] == pytest.approx(1000.0)   # clamped above
