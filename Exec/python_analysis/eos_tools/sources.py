"""Raw-source manifest handling: checksum verification and scripted fetch.

The manifest (`Exec/testing/EOS-Table/data/raw/sources.yaml`) records every
raw building-block dataset: citation, origin URL/DOI, retrieval date, sha256
per file, license note and acquisition mode (scripted | manual | digitized).
Table provenance headers point back at it.

`verify` recomputes every checksum (the SS1 gate; always runnable).
`fetch` downloads only the scriptable subset, for files not already present,
then verifies them. Manual/digitized entries are reported, never fetched.
"""

import hashlib
import os
import sys
import urllib.request

import yaml


def repo_root():
    # this file lives at Exec/python_analysis/eos_tools/sources.py
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "..", ".."))


def default_manifest():
    return os.path.join(repo_root(), "Exec", "testing", "EOS-Table", "data",
                        "raw", "sources.yaml")


def sha256_of(path, bufsize=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_manifest(path):
    with open(path) as f:
        doc = yaml.safe_load(f)
    if not isinstance(doc, dict) or "sources" not in doc:
        raise ValueError("%s: expected a top-level 'sources' list" % path)
    return doc


def verify(manifest_path, quiet=False):
    """Check every file in the manifest. Returns the number of failures."""
    base = os.path.dirname(os.path.abspath(manifest_path))
    doc = load_manifest(manifest_path)
    n_fail = 0
    for src in doc["sources"]:
        for fent in src.get("files", []):
            path = os.path.join(base, fent["path"])
            if not os.path.exists(path):
                print("MISSING   %-12s %s" % (src["id"], fent["path"]))
                n_fail += 1
                continue
            got = sha256_of(path)
            if got != fent["sha256"]:
                print("MISMATCH  %-12s %s\n  expected %s\n  got      %s"
                      % (src["id"], fent["path"], fent["sha256"], got))
                n_fail += 1
            elif not quiet:
                print("OK        %-12s %s" % (src["id"], fent["path"]))
    print("verify: %d failure(s)" % n_fail)
    return n_fail


def fetch(manifest_path):
    """Download scriptable sources whose files are absent, then verify them.

    Returns the number of failures (download errors or checksum mismatches).
    """
    base = os.path.dirname(os.path.abspath(manifest_path))
    doc = load_manifest(manifest_path)
    n_fail = 0
    for src in doc["sources"]:
        if src.get("acquisition") != "scripted":
            print("skip (%s) %-12s — acquire per its notes/citation"
                  % (src.get("acquisition", "?"), src["id"]))
            continue
        for fent in src.get("files", []):
            path = os.path.join(base, fent["path"])
            if os.path.exists(path):
                continue
            url = fent.get("url") or src.get("url")
            if not url:
                print("NO-URL    %-12s %s" % (src["id"], fent["path"]))
                n_fail += 1
                continue
            print("fetching  %-12s %s\n  <- %s" % (src["id"], fent["path"], url))
            os.makedirs(os.path.dirname(path), exist_ok=True)
            tmp = path + ".part"
            try:
                with urllib.request.urlopen(url, timeout=120) as r, \
                        open(tmp, "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
                os.replace(tmp, path)
            except Exception as ex:
                print("FETCH-FAIL %-12s %s: %s" % (src["id"], fent["path"], ex))
                if os.path.exists(tmp):
                    os.remove(tmp)
                n_fail += 1
                continue
            got = sha256_of(path)
            if got != fent["sha256"]:
                print("MISMATCH  %-12s %s (post-fetch)\n  expected %s\n  got      %s"
                      % (src["id"], fent["path"], fent["sha256"], got))
                n_fail += 1
    print("fetch: %d failure(s)" % n_fail)
    return n_fail


def main_sources(args):
    manifest = args.manifest or default_manifest()
    if not os.path.exists(manifest):
        print("manifest not found: %s" % manifest)
        sys.exit(2)
    n_fail = 0
    if args.fetch:
        n_fail += fetch(manifest)
    if args.verify or not args.fetch:
        n_fail += verify(manifest, quiet=args.quiet)
    sys.exit(1 if n_fail else 0)
