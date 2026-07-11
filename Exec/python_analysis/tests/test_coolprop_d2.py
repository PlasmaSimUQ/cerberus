"""CoolProp/Richardson-2014 fluid-D2 wrapper (skipped if CoolProp absent).

The reference-state closure numbers here ARE two of the SS2 exit gates
(doc/eos_creation_plan.md SS2).
"""

import numpy as np
import pytest

CP = pytest.importorskip("CoolProp.CoolProp")

from eos_tools.models.coolprop_d2 import CoolPropD2  # noqa: E402


def test_reference_state_density():
    # SS2 gate: rho(20 K, 1 bar) = 0.171 g/cc +- 0.5%
    rho = CP.PropsSI("D", "T", 20.0, "P", 1.0e5, "Deuterium") * 1e-3
    assert abs(rho / 0.171 - 1.0) < 5e-3, rho


def test_reference_state_sound_speed():
    # SS2 gate: liquid cs within 10% of ~1.1 km/s
    cs = CP.PropsSI("A", "T", 20.0, "P", 1.0e5, "Deuterium")
    assert 0.9 * 1100.0 < cs < 1.1 * 1100.0, cs


def test_wrapper_masks_and_consistency():
    m = CoolPropD2()
    # exact 1-bar liquid density: rho = 0.171 itself sits just inside the
    # two-phase dome at 20 K (rho_satL ~ 171.8 kg/m3) and must be masked
    rho_1bar = CP.PropsSI("D", "T", 20.0, "P", 1.0e5, "Deuterium") * 1e-3
    rho = np.array([rho_1bar, 0.171, 0.05])
    T = np.array([20.0, 5.0, 300.0])  # 5 K is below Tmin -> invalid
    out = m.eval(rho, T)
    assert bool(out["valid"][0]) and not bool(out["valid"][1])
    # a = e - T s by construction
    i = out["valid"]
    assert np.allclose(out["a"][i], out["e"][i] - T[i] * out["s"][i],
                       rtol=1e-12)


def test_wrapper_pressure_positive_and_cgs():
    m = CoolPropD2()
    rho_1bar = CP.PropsSI("D", "T", 20.0, "P", 1.0e5, "Deuterium") * 1e-3
    out = m.eval(rho_1bar, 20.0)
    # ~1 bar in barye at the reference point (loose: just the unit chain)
    assert 1e5 < out["p"] < 1e8
    assert out["cs"] > 5e4  # cm/s
