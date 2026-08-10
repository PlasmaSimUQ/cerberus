"""Common-energy-reference plumbing (HANDOFF addendum): sample_bilinear
must reproduce the C++ reader's evaluation convention (bilinear in the
log10 axes, linear in the block value), since --e-ref-state anchors the
mixture-set gauge to what the solver will actually read."""

import numpy as np
import pytest

from eos_tools.condition import sample_bilinear


def _grid():
    lrho = np.linspace(-4.0, 1.0, 11)
    lT = np.linspace(2.0, 6.0, 9)
    return lrho, lT


def test_bilinear_exact_on_plane():
    # a function linear in (lrho, lT) is reproduced exactly everywhere
    lrho, lT = _grid()
    F = 3.0 * lrho[:, None] - 2.0 * lT[None, :] + 7.0
    for rho, T in [(1e-3, 1e4), (0.5, 3.7e3), (10.0 ** -3.95, 10.0 ** 5.9)]:
        want = 3.0 * np.log10(rho) - 2.0 * np.log10(T) + 7.0
        assert sample_bilinear(lrho, lT, F, rho, T) == pytest.approx(
            want, rel=1e-13)


def test_node_and_endpoint_values():
    lrho, lT = _grid()
    rng = np.random.default_rng(7)
    F = rng.normal(size=(len(lrho), len(lT)))
    # exact at interior nodes and at the top corners (i,j clamped to n-2)
    assert sample_bilinear(lrho, lT, F, 10.0 ** lrho[3], 10.0 ** lT[5]) \
        == pytest.approx(F[3, 5], rel=1e-12)
    assert sample_bilinear(lrho, lT, F, 10.0 ** lrho[-1], 10.0 ** lT[-1]) \
        == pytest.approx(F[-1, -1], rel=1e-12)


def test_linear_in_value_not_loglog():
    # midpoint of a cell: mean of the four corner VALUES (the reader
    # interpolates values linearly; a log-log sampler would return the
    # geometric mean and fail this)
    lrho, lT = _grid()
    F = np.ones((len(lrho), len(lT)))
    F[4, 4] = 100.0
    x = 10.0 ** (0.5 * (lrho[4] + lrho[5]))
    y = 10.0 ** (0.5 * (lT[4] + lT[5]))
    assert sample_bilinear(lrho, lT, F, x, y) == pytest.approx(
        (100.0 + 1.0 + 1.0 + 1.0) / 4.0, rel=1e-12)


def test_out_of_range_raises():
    lrho, lT = _grid()
    F = np.zeros((len(lrho), len(lT)))
    with pytest.raises(ValueError, match="outside the table axes"):
        sample_bilinear(lrho, lT, F, 1e-5, 1e4)  # below rho axis
    with pytest.raises(ValueError, match="outside the table axes"):
        sample_bilinear(lrho, lT, F, 1.0, 10.0 ** 6.1)  # above T axis
