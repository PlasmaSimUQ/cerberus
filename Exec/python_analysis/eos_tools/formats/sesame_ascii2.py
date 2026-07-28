"""SESAME ASCII2 ingest (LANL sesame-unc distribution).

File format: LA-UR-25-26256 ("SESAME: ASCII2 File Format", 2025); table
payload structures: LA-UR-25-26258 ("SESAME ... with Extensions for
Multi-phase Representations", 2025). Facts this reader relies on:

- Line 1: ``Version 2.0`` (case-insensitive label).
- Table header: one line of 7 whitespace-delimited positive integers,
  ``file_no mat_id table_id n_words create_date update_date version``
  (file_no 0 starts a material, 1 continues it; no end-of-file record).
- Comment tables (table_id 101-199, 10101-10199): n_words = character
  count excluding end-of-line characters; 80 chars/line except the last.
- Numeric tables: n_words whitespace-delimited floats, <= 5 per line, any
  C-parseable float syntax. The reader tokenizes on whitespace and counts
  to n_words, ignoring line structure entirely.
- 2-D EOS payload (301/303/304/305/311 and phase-specific kin):
  ``NR NT rho[NR] T[NT] P[NR*NT] U[NR*NT] (A[NR*NT])`` with density
  varying FASTEST (Fortran/column-major: value(i,j) at index i + NR*j).
  The Helmholtz array A is absent from most legacy tables — detected from
  n_words, never assumed (LA-UR-25-26258 footnote 11).
- Units: rho Mg/m^3 (== g/cc), T K, P GPa, U and A MJ/kg.

Streaming: the production sesame-unc.ascii2 is ~244 MB; everything here
does a single sequential pass and materialises only the requested tables.
"""

import os

import numpy as np

from ..constants import GPA_CGS

MJKG_CGS = 1.0e10  # MJ/kg -> erg/g

# 2-D (rho, T) EOS families this reader can hand to the eostab pipeline
TWO_D_TABLES = (301, 303, 304, 305, 311)


def _is_comment(tid):
    return 101 <= tid <= 199 or 10101 <= tid <= 10199


def _read_comment(f, n_chars):
    """Read comment lines until n_chars non-EOL characters are consumed."""
    lines, seen = [], 0
    while seen < n_chars:
        line = f.readline()
        if not line:
            raise ValueError(
                "EOF inside comment table (%d of %d chars read)" % (seen, n_chars))
        s = line.rstrip("\r\n")
        lines.append(s)
        seen += len(s)
    if seen != n_chars:
        raise ValueError(
            "comment table overran its header count (%d != %d); file "
            "corrupt or a non-80-column comment line" % (seen, n_chars))
    return lines


def _read_floats(f, n_words):
    """Read exactly n_words whitespace-delimited floats."""
    vals = []
    while len(vals) < n_words:
        line = f.readline()
        if not line:
            raise ValueError(
                "EOF inside numeric table (%d of %d words read)" % (len(vals), n_words))
        vals.extend(float(t) for t in line.split())
    if len(vals) != n_words:
        raise ValueError(
            "numeric table overran its header count (%d != %d)" % (len(vals), n_words))
    return np.array(vals)


def _skip_floats(f, n_words):
    """Skip a numeric payload without float conversion (cheap pass)."""
    seen = 0
    while seen < n_words:
        line = f.readline()
        if not line:
            raise ValueError(
                "EOF inside numeric table (%d of %d words skipped)" % (seen, n_words))
        seen += len(line.split())
    if seen != n_words:
        raise ValueError(
            "numeric table overran its header count (%d != %d)" % (seen, n_words))


def _parse_header(line):
    """Header line -> (file_no, mat_id, table_id, n_words) or None."""
    t = line.split()
    if len(t) != 7 or t[0] not in ("0", "1"):
        return None
    if not all(x.isdigit() for x in t):
        return None
    return int(t[0]), int(t[1]), int(t[2]), int(t[3])


def _open_checked(path):
    f = open(path)
    first = f.readline()
    t = first.split()
    if len(t) != 2 or t[0].lower() != "version":
        f.close()
        raise ValueError("%s: not a SESAME ASCII2 file (first line %r)"
                         % (path, first.strip()))
    return f


def _walk(path):
    """Yield (f, mat_id, table_id, n_words) for every table, payload unread.

    The consumer MUST consume or skip the payload before advancing the
    generator (position-based parsing: after a payload the next line is by
    construction the next header).
    """
    f = _open_checked(path)
    try:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.strip():
                continue
            h = _parse_header(line)
            if h is None:
                raise ValueError("%s: expected a table header, got %r"
                                 % (path, line.rstrip()))
            yield f, h[1], h[2], h[3]
    finally:
        f.close()


def index(path):
    """One cheap pass: {mat_id: {'name': <101 material field>, 'tables':
    {table_id: n_words}}}, insertion-ordered as found in the file."""
    mats = {}
    for f, mat, tid, nw in _walk(path):
        m = mats.setdefault(mat, {"name": None, "tables": {}})
        m["tables"][tid] = nw
        if _is_comment(tid):
            lines = _read_comment(f, nw)
            if tid == 101 and m["name"] is None:
                text = " ".join(lines)
                key = "material."
                i = text.lower().find(key)
                if i >= 0:
                    j = text.find("/", i)
                    m["name"] = text[i + len(key):j if j > 0 else None].strip()
        else:
            _skip_floats(f, nw)
    return mats


def read_sesame(path, mat_id, table_id):
    """Extract one 2-D EOS table for one material, in SESAME units.

    Returns a dict:
      rho[NR] (Mg/m^3 == g/cc), T[NT] (K),
      p[NR,NT] (GPa), e[NR,NT] (MJ/kg), a[NR,NT] (MJ/kg) or None,
      zbar, abar, rho0, bs0 (from the 201 table; None if absent),
      comment101 (str), has_helmholtz (bool),
      mat_id, table_id, src (basename).

    Arrays are returned TRANSPOSED to (i_rho, j_T) C-order — i.e. p[i, j]
    at (rho[i], T[j]) — matching the eos_tools (Nr, Nt) convention; the
    SESAME stream stores rho-fastest and the reshape accounts for it.
    """
    if table_id not in TWO_D_TABLES:
        raise ValueError(
            "table %d is not a 2-D (rho,T) EOS table; supported: %s "
            "(306 is a 1-D cold curve, 4xx are phase boundaries — neither "
            "can populate a (rho,T) .eostab)" % (table_id, TWO_D_TABLES))
    out = {"mat_id": mat_id, "table_id": table_id,
           "src": os.path.basename(path), "zbar": None, "abar": None,
           "rho0": None, "bs0": None, "comment101": "", "a": None}
    found = False
    for f, mat, tid, nw in _walk(path):
        if mat != mat_id:
            if _is_comment(tid):
                _read_comment(f, nw)
            else:
                _skip_floats(f, nw)
            continue
        if tid == 101:
            out["comment101"] = " ".join(_read_comment(f, nw))
        elif _is_comment(tid):
            _read_comment(f, nw)
        elif tid == 201:
            v = _read_floats(f, nw)
            out["zbar"], out["abar"] = v[0], v[1]
            out["rho0"] = v[2]
            if nw > 3:
                out["bs0"] = v[3]
        elif tid == table_id:
            v = _read_floats(f, nw)
            nr, nt = int(v[0]), int(v[1])
            if nr < 2 or nt < 2:
                raise ValueError("material %d table %d: degenerate grid %dx%d"
                                 % (mat_id, table_id, nr, nt))
            nfun = (nw - 2 - nr - nt) / (nr * nt)
            if nfun not in (2.0, 3.0):
                raise ValueError(
                    "material %d table %d: n_words=%d inconsistent with "
                    "%dx%d grid (implies %.3f function arrays; expected 2 "
                    "[P,U] or 3 [P,U,A])" % (mat_id, table_id, nw, nr, nt, nfun))
            nfun = int(nfun)
            out["rho"] = v[2:2 + nr]
            out["T"] = v[2 + nr:2 + nr + nt]
            k = 2 + nr + nt
            # SESAME stores rho-fastest: reshape (nt, nr) then transpose
            arrs = [v[k + a * nr * nt:k + (a + 1) * nr * nt]
                    .reshape(nt, nr).T.copy() for a in range(nfun)]
            out["p"], out["e"] = arrs[0], arrs[1]
            if nfun == 3:
                out["a"] = arrs[2]
            out["has_helmholtz"] = (nfun == 3)
            found = True
        else:
            _skip_floats(f, nw)
    if not found:
        raise ValueError("material %d has no table %d in %s"
                         % (mat_id, table_id, path))
    for ax in ("rho", "T"):
        d = np.diff(out[ax])
        if not np.all(d >= 0.0):
            raise ValueError("material %d table %d: %s axis not "
                             "monotonically increasing" % (mat_id, table_id, ax))
    return out


def sesame_to_cgs(raw):
    """SESAME units -> the standard scattered-points dict (CGS).

    rho Mg/m^3 == g/cc (x1); P GPa -> erg/cc (x1e10); U MJ/kg -> erg/g
    (x1e10). The rectangular grid is flattened to scattered points
    (grids.regrid groups them back into isochores by unique rho).
    """
    rho, T = raw["rho"], raw["T"]
    R = np.repeat(rho, len(T))
    Tt = np.tile(T, len(rho))
    return {
        "rho": R,
        "T": Tt,
        "p": raw["p"].ravel(order="C") * GPA_CGS,
        "e": raw["e"].ravel(order="C") * MJKG_CGS,
    }
