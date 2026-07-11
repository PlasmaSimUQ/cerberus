"""Exact-blend identities (models/composite.py) + alignment recovery."""

import numpy as np

from eos_tools.models.base import HelmholtzModel
from eos_tools.models.composite import (LogRamp, align_constants, blend,
                                        shifted)


class PowerLawModel(HelmholtzModel):
    """a = C rho^alpha T^beta — closed-form p, e, s for reference."""

    def __init__(self, C, alpha, beta):
        self.C, self.al, self.be = C, alpha, beta

    def a(self, rho, T):
        return self.C * np.asarray(rho, float) ** self.al \
            * np.asarray(T, float) ** self.be

    def out(self, rho, T):
        rho, T = np.asarray(rho, float), np.asarray(T, float)
        a = self.a(rho, T)
        return {"a": a, "p": rho * self.al * a,
                "s": -self.be * a / T, "e": a * (1.0 - self.be)}


def test_blend_of_identical_sources_is_identity():
    m = PowerLawModel(3.0e8, 0.5, 1.3)
    rho = np.geomspace(1e-3, 1.0, 12)[:, None]
    T = np.geomspace(30.0, 3000.0, 10)[None, :]
    o = m.out(rho * np.ones_like(T), T * np.ones_like(rho))
    ramp = LogRamp(0.0, 0.1)
    res = blend(o, o, ramp.w(rho), dw_drho=ramp.dw_dval(rho), rho=rho)
    for k in ("a", "p", "e", "s"):
        assert np.allclose(res[k], o[k], rtol=1e-14), k


def test_blend_matches_differentiated_helmholtz():
    """Exact formulas == numerically differentiating the blended a."""
    m1 = PowerLawModel(3.0e8, 0.5, 1.3)
    m2 = PowerLawModel(1.0e8, 1.1, 0.7)
    rampR = LogRamp(-1.0, 0.15)   # in rho
    rampT = LogRamp(2.3, 0.12)    # in T

    class Blended(HelmholtzModel):
        def a(self, rho, T):
            w = rampR.w(rho) * rampT.w(T)
            return w * m1.a(rho, T) + (1.0 - w) * m2.a(rho, T)

    rho = np.array([0.05, 0.1, 0.2])
    T = np.array([150.0, 200.0, 300.0])
    o1, o2 = m1.out(rho, T), m2.out(rho, T)
    w = rampR.w(rho) * rampT.w(T)
    dw_drho = rampR.dw_dval(rho) * rampT.w(T)
    dw_dT = rampR.w(rho) * rampT.dw_dval(T)
    res = blend(o1, o2, w, dw_drho=dw_drho, dw_dT=dw_dT, rho=rho, T=T)

    ref = Blended()
    assert np.allclose(res["p"], ref.p(rho, T), rtol=1e-6)
    assert np.allclose(res["s"], ref.s(rho, T), rtol=1e-6)
    assert np.allclose(res["e"], ref.e(rho, T), rtol=1e-6)


def test_align_recovers_planted_offsets():
    rng = np.random.default_rng(7)
    e_ref = rng.uniform(1e9, 2e9, 40)
    s_ref = rng.uniform(1e7, 2e7, 40)
    cE_true, cS_true = 3.3e8, -4.4e6
    cE, cS, std_e, std_s = align_constants(e_ref, s_ref,
                                           e_ref - cE_true, s_ref - cS_true)
    assert abs(cE / cE_true - 1.0) < 1e-12
    assert abs(cS / cS_true - 1.0) < 1e-12
    assert std_e < 1e-3 and std_s < 1e-9


def test_shift_leaves_p_and_moves_e_s():
    m = PowerLawModel(2.0e8, 0.8, 1.1)
    rho, T = np.array([0.1]), np.array([250.0])
    o = m.out(rho, T)
    sh = shifted(o, 1e9, 1e6, T)
    assert np.allclose(sh["p"], o["p"])
    assert np.allclose(sh["e"], o["e"] + 1e9)
    assert np.allclose(sh["s"], o["s"] + 1e6)
    assert np.allclose(sh["a"], o["a"] + 1e9 - T * 1e6)
