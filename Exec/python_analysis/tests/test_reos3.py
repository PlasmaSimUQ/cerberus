"""H-REOS.3 reader invariants + isotope scaling (uses the committed file)."""

import os

import numpy as np
import pytest

from eos_tools.constants import GPA_CGS, M_D, M_H
from eos_tools.formats.reos3 import read_reos3, reos3_to_deuterium

PATH = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "testing", "EOS-Table", "data",
    "raw", "reos3", "table2_HREOS3.dat"))


@pytest.fixture(scope="module")
def raw():
    if not os.path.exists(PATH):
        pytest.skip("REOS.3 raw table not present")
    return read_reos3(PATH)


def test_shape_and_ranges(raw):
    assert len(raw["rho"]) == 4410
    assert raw["n_isotherms"] == 42
    assert raw["T"].min() == 60.0 and raw["T"].max() == 1.0e7
    assert np.all(raw["p"] > 0.0)
    # CGS: 1 GPa = 1e10 erg/cc round trip on a spot value
    assert raw["p"].max() > 1e5 * GPA_CGS  # multi-TPa top end


def test_deuterium_scaling(raw):
    pts = reos3_to_deuterium(raw)
    s = M_D / M_H
    assert np.allclose(pts["rho"], raw["rho"] * s)
    assert np.allclose(pts["e"], raw["e"] / s)
    assert np.allclose(pts["p"], raw["p"])
    # equal nuclear number density preserved: rho/m per basis identical
    assert np.allclose(pts["rho"] / M_D, raw["rho"] / M_H)


def test_reader_rejects_corruption(tmp_path, raw):
    import numpy as np
    bad = tmp_path / "bad.dat"
    rho, T = raw["rho"], raw["T"]
    P = (raw["p"] / GPA_CGS).copy()
    u = raw["e"] / 1e10
    P[5] = -1.0
    np.savetxt(bad, np.column_stack([rho, T, P, u]), fmt="%.6e")
    with pytest.raises(ValueError):
        read_reos3(str(bad))
