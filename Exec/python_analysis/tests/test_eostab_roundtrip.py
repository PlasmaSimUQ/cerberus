"""write_eostab -> read_eostab must round-trip axes, blocks and provenance."""

import numpy as np

from eos_tools.formats.eostab import read_eostab, write_eostab


def make_blocks(nr, nt, ne, npp):
    rng = np.random.default_rng(42)
    blocks = {name: rng.uniform(1.0, 2.0, (nr, nt))
              for name in ("p", "e", "dpdT", "dpdrho", "cv", "dedrho")}
    blocks["hull"] = (rng.uniform(0, 1, (nr, nt)) > 0.3).astype(float)
    blocks["T_of_e"] = rng.uniform(1e3, 1e7, (nr, ne))
    blocks["T_of_p"] = rng.uniform(1e3, 1e7, (nr, npp))
    return blocks


def test_roundtrip(tmp_path):
    nr, nt, ne, npp = 7, 5, 6, 4
    lrho = np.linspace(-3.0, 1.0, nr)
    lT = np.linspace(3.0, 7.0, nt)
    le = np.linspace(10.0, 14.0, ne)
    lp = np.linspace(6.0, 16.0, npp)
    blocks = make_blocks(nr, nt, ne, npp)
    prov = [("material", "test"), ("units", "cgs"), ("e_shift", "0.0")]
    path = str(tmp_path / "t.eostab")
    write_eostab(path, prov, lrho, lT, blocks, {"le": le, "lp": lp})

    prov2, meta, blocks2 = read_eostab(path)
    assert prov2["material"] == "test"
    assert np.allclose(meta["lrho"], lrho)
    assert np.allclose(meta["lT"], lT)
    assert np.allclose(meta["le"], le)
    assert np.allclose(meta["lp"], lp)
    for name, arr in blocks.items():
        assert blocks2[name].shape == arr.shape, name
        # writer uses %.10e -> relative round-off ~1e-10
        assert np.allclose(blocks2[name], arr, rtol=1e-9, atol=0.0), name


def test_writer_rejects_nonfinite(tmp_path):
    nr, nt = 4, 4
    lrho = np.linspace(-1, 1, nr)
    lT = np.linspace(3, 5, nt)
    blocks = make_blocks(nr, nt, nt, nt)
    blocks["p"][2, 2] = np.nan
    try:
        write_eostab(str(tmp_path / "bad.eostab"),
                     [("material", "x")], lrho, lT, blocks,
                     {"le": lT, "lp": lT})
    except AssertionError as ex:
        assert "non-finite" in str(ex)
    else:
        raise AssertionError("writer accepted a NaN block")
