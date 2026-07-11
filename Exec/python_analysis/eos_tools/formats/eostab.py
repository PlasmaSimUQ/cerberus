"""The canonical ``.eostab`` writer/reader.

.eostab format spec v1 — FROZEN 2026-07-07
(mirror of Exec/testing/EOS-Table/README.md; do not change one without the
other):

    EOSTAB 1                          # magic + format version
    # provenance: free-form 'key: value' lines, '#' comments allowed
    material:    D
    source:      <origin, retrieval date>
    generator:   eos_table_prep.py <git-sha>, run <date>
    composition: A=2.014 Z=1
    units:       cgs                  # values stored DIMENSIONAL (CGS):
                                      #   rho g/cc, T K, p erg/cc, e erg/g
                                      # Cerberus nondimensionalises at load
    e_shift:     <erg/g>              # constant added to raw specific
                                      # internal energy so min(e) > 0 on the
                                      # hull (zero-point is arbitrary; must
                                      # be used consistently, and is)
    conditioning: cv_floor=<v> monotonised=<n> cv_floored=<n> maxwell=<status>
    # grid (uniform in log10); axes reconstructed by reader, never stored
    grid: n_rho=<Nr> n_T=<Nt>
    lrho: <log10 rho_min> <log10 rho_max>
    lT:   <log10 T_min>   <log10 T_max>
    # inverse-map axes (required iff T_of_e / T_of_p blocks present)
    le:   <log10 e_min> <log10 e_max>  n_e=<Ne>
    lp:   <log10 p_min> <log10 p_max>  n_p=<Np>
    # data blocks: 'block: <name>' then the values, whitespace separated.
    # storage order: i_rho outer, j inner  (idx = i*n_T + j; inverse maps
    # idx = i*n_e + j / i*n_p + j)
    block: p          # pressure
    block: e          # specific internal energy (shifted)
    block: dpdT       # (dP/dT)|rho     smoothed  -- outputs only
    block: dpdrho     # (dP/drho)|T     smoothed/monotonised
    block: cv         # (de/dT)|rho     smoothed, floored > 0
    block: dedrho     # (de/drho)|T     smoothed
    block: T_of_e     # T(rho, e) inverse map  -- Newton seeds (plan D5/D8)
    block: T_of_p     # T(rho, p) inverse map
    block: hull       # REQUIRED: 1.0 inside source hull, 0.0 filled cell
"""

import numpy as np


def write_eostab(path, prov, lrho, lT, blocks, inv_axes):
    """Write a .eostab per the frozen spec (see module docstring).

    prov: list of (key, value) provenance lines, in order.
    lrho/lT: 1-D log10 axes (uniform).
    blocks: dict name -> 2-D array; forward blocks shaped (Nr, Nt),
            T_of_e (Nr, Ne), T_of_p (Nr, Np).
    inv_axes: dict with 'le' (1-D log10 e axis), 'lp' (1-D log10 p axis).
    """
    nr, nt = len(lrho), len(lT)
    order = ["p", "e", "dpdT", "dpdrho", "cv", "dedrho", "T_of_e", "T_of_p", "hull"]
    with open(path, "w") as f:
        f.write("EOSTAB 1\n")
        for k, v in prov:
            f.write("%s: %s\n" % (k, v))
        f.write("grid: n_rho=%d n_T=%d\n" % (nr, nt))
        f.write("lrho: %.12g %.12g\n" % (lrho[0], lrho[-1]))
        f.write("lT: %.12g %.12g\n" % (lT[0], lT[-1]))
        le, lp = inv_axes["le"], inv_axes["lp"]
        f.write("le: %.12g %.12g n_e=%d\n" % (le[0], le[-1], len(le)))
        f.write("lp: %.12g %.12g n_p=%d\n" % (lp[0], lp[-1], len(lp)))
        for name in order:
            arr = np.asarray(blocks[name])
            assert arr.shape[0] == nr, (name, arr.shape)
            assert np.all(np.isfinite(arr)), "non-finite values in block " + name
            f.write("block: %s\n" % name)
            flat = arr.ravel(order="C")  # idx = i*ncol + j
            for i in range(0, flat.size, 6):
                f.write(" ".join("%.10e" % v for v in flat[i:i + 6]) + "\n")
    print("wrote %s  (%d x %d, %d blocks)" % (path, nr, nt, len(order)))


def read_eostab(path):
    """Minimal reader (QA use; prototype of the C++ parser)."""
    prov, blocks = {}, {}
    with open(path) as f:
        tok = f.read().split("\n")
    assert tok[0].startswith("EOSTAB 1"), "bad magic"
    i, cur = 1, None
    vals = []
    meta = {}
    while i < len(tok):
        line = tok[i].strip()
        i += 1
        if not line or line.startswith("#"):
            continue
        if line.startswith("block:"):
            if cur:
                blocks[cur] = np.array(vals)
            cur, vals = line.split()[1], []
        elif cur is not None:
            vals.extend(float(x) for x in line.split())
        else:
            k, v = line.split(":", 1)
            prov[k.strip()] = v.strip()
    if cur:
        blocks[cur] = np.array(vals)
    g = dict(p.split("=") for p in prov["grid"].split())
    nr, nt = int(g["n_rho"]), int(g["n_T"])
    a, b = (float(x) for x in prov["lrho"].split())
    meta["lrho"] = np.linspace(a, b, nr)
    a, b = (float(x) for x in prov["lT"].split())
    meta["lT"] = np.linspace(a, b, nt)
    for ax, nkey in (("le", "n_e"), ("lp", "n_p")):
        parts = prov[ax].split()
        n = int(dict(p.split("=") for p in parts if "=" in p)[nkey])
        meta[ax] = np.linspace(float(parts[0]), float(parts[1]), n)
    for name in blocks:
        ncol = {"T_of_e": len(meta["le"]), "T_of_p": len(meta["lp"])}.get(name, nt)
        blocks[name] = blocks[name].reshape(nr, ncol)
    return prov, meta, blocks
