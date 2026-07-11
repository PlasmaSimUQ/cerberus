"""Manifest verify/fetch logic (no network: fetch only exercises skip paths)."""

import os

import yaml

from eos_tools import sources


def write_manifest(tmp_path, entries):
    man = tmp_path / "sources.yaml"
    with open(man, "w") as f:
        yaml.safe_dump({"sources": entries}, f)
    return str(man)


def test_verify_ok_and_mismatch(tmp_path):
    data = tmp_path / "blob.txt"
    data.write_text("hello eos\n")
    good = sources.sha256_of(str(data))

    man = write_manifest(tmp_path, [{
        "id": "blob", "acquisition": "manual",
        "files": [{"path": "blob.txt", "sha256": good}],
    }])
    assert sources.verify(man) == 0

    data.write_text("tampered\n")
    assert sources.verify(man) == 1


def test_verify_missing_file(tmp_path):
    man = write_manifest(tmp_path, [{
        "id": "ghost", "acquisition": "scripted",
        "files": [{"path": "not_there.bin", "sha256": "0" * 64}],
    }])
    assert sources.verify(man) == 1


def test_fetch_skips_manual_and_present(tmp_path):
    data = tmp_path / "have.txt"
    data.write_text("present\n")
    man = write_manifest(tmp_path, [
        {"id": "manual-src", "acquisition": "digitized",
         "files": [{"path": "x.csv", "sha256": "0" * 64}]},
        {"id": "have-src", "acquisition": "scripted",
         "files": [{"path": "have.txt",
                    "sha256": sources.sha256_of(str(data))}]},
    ])
    # digitized entry skipped, present file not re-fetched -> no failures
    assert sources.fetch(man) == 0


def test_fetch_reports_missing_url(tmp_path):
    man = write_manifest(tmp_path, [{
        "id": "nourl", "acquisition": "scripted",
        "files": [{"path": "y.bin", "sha256": "0" * 64}],
    }])
    assert sources.fetch(man) == 1


def test_default_manifest_path_inside_repo():
    p = sources.default_manifest()
    assert p.endswith(os.path.join("EOS-Table", "data", "raw", "sources.yaml"))
