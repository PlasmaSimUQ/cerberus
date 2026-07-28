"""SESAME ASCII2 parser invariants on a synthetic in-test fixture.

The fixture exercises the format facts the reader relies on (LA-UR-25-26256):
free whitespace (spaces + tabs), any float syntax, 5-per-line data with a
ragged last line, comment tables counted in non-EOL characters, optional
Helmholtz array, rho-fastest storage order.
"""

import numpy as np
import pytest

from eos_tools.constants import GPA_CGS
from eos_tools.formats.sesame_ascii2 import (
    MJKG_CGS, index, read_sesame, sesame_to_cgs)

RHO = np.array([1.0, 2.0, 3.0, 4.0])       # Mg/m^3
TEMP = np.array([10.0, 20.0, 30.0])        # K
NR, NT = len(RHO), len(TEMP)


def p_of(i, j):
    return RHO[i] + 100.0 * TEMP[j]        # GPa


def e_of(i, j):
    return 7.0 * RHO[i] + TEMP[j]          # MJ/kg


def a_of(i, j):
    return e_of(i, j) - 0.5 * TEMP[j]      # MJ/kg


def stream_301(nfun):
    """Payload in SESAME storage order: rho varies fastest (i + NR*j)."""
    vals = [NR, NT] + list(RHO) + list(TEMP)
    funcs = [p_of, e_of, a_of][:nfun]
    for fn in funcs:
        for j in range(NT):
            for i in range(NR):
                vals.append(fn(i, j))
    return vals


def five_per_line(vals, sep=" "):
    lines = []
    for k in range(0, len(vals), 5):
        lines.append(sep.join("%.17g" % v for v in vals[k:k + 5]))
    return "\n".join(lines)


@pytest.fixture(scope="module")
def fixture(tmp_path_factory):
    c1a = "material. testium (z=1.0, a=2.0) /source. pytest /comp. X /".ljust(80)
    c1b = "second line"                   # only the last line may be short
    n1 = len(c1a) + len(c1b)              # EOL characters not counted
    generic = "a generic comment line under 80 chars"
    body = ["Version 2.0"]
    # material 42: 101, 201, 301 with Helmholtz; tabs as separators
    body.append(" 0 42 101 %d 20240101 20240101 1" % n1)
    body.append(c1a)
    body.append(c1b)
    body.append(" 1 42 201 5 20240101 20240101 1")
    body.append(five_per_line([1.0, 2.0, 3.5, 110.0, 0.0]))
    v3 = stream_301(3)
    body.append(" 1\t42\t301 %d 20240101 20240101 1" % len(v3))
    body.append(five_per_line(v3, sep="\t"))
    # material 43: generic comment + 301 without Helmholtz, terse floats
    body.append(" 0 43 102 %d 20240101 20240101 1" % len(generic))
    body.append(generic)
    v2 = stream_301(2)
    body.append(" 1 43 301 %d 20240101 20240101 1" % len(v2))
    body.append(five_per_line(v2))
    path = tmp_path_factory.mktemp("ses") / "mini.ascii2"
    path.write_text("\n".join(body) + "\n")
    return str(path)


def test_index(fixture):
    mats = index(fixture)
    assert set(mats) == {42, 43}
    assert mats[42]["name"] == "testium (z=1.0, a=2.0)"
    assert set(mats[42]["tables"]) == {101, 201, 301}
    assert mats[43]["name"] is None       # no 101 table
    assert set(mats[43]["tables"]) == {102, 301}


def test_read_with_helmholtz_and_transpose(fixture):
    raw = read_sesame(fixture, 42, 301)
    assert raw["has_helmholtz"]
    assert raw["zbar"] == 1.0 and raw["abar"] == 2.0
    assert raw["rho0"] == 3.5 and raw["bs0"] == 110.0
    assert "testium" in raw["comment101"]
    assert np.array_equal(raw["rho"], RHO)
    assert np.array_equal(raw["T"], TEMP)
    for i in range(NR):
        for j in range(NT):
            assert raw["p"][i, j] == p_of(i, j)     # transpose is exact
            assert raw["e"][i, j] == e_of(i, j)
            assert raw["a"][i, j] == a_of(i, j)


def test_read_two_function_table(fixture):
    raw = read_sesame(fixture, 43, 301)
    assert not raw["has_helmholtz"]
    assert raw["a"] is None
    assert raw["zbar"] is None            # no 201 table for this material
    assert raw["p"][2, 1] == p_of(2, 1)


def test_cgs_conversion(fixture):
    raw = read_sesame(fixture, 42, 301)
    pts = sesame_to_cgs(raw)
    assert len(pts["rho"]) == NR * NT
    k = np.flatnonzero((pts["rho"] == RHO[1]) & (pts["T"] == TEMP[2]))[0]
    assert pts["p"][k] == p_of(1, 2) * GPA_CGS
    assert pts["e"][k] == e_of(1, 2) * MJKG_CGS


def test_errors(fixture):
    with pytest.raises(ValueError, match="no table 304"):
        read_sesame(fixture, 42, 304)
    with pytest.raises(ValueError, match="not a 2-D"):
        read_sesame(fixture, 42, 306)
    with pytest.raises(ValueError, match="not a SESAME ASCII2"):
        index(__file__)
