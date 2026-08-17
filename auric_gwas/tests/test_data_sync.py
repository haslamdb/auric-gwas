"""Guard the three vendored family-definition files against silent drift from staph_snp_db.

The manifest check (test_vendored_files_match_manifest) is self-contained and always runs.
The source check (test_source_matches_manifest) skips unless a staph_snp_db checkout is present.

Run: python -m pytest auric_gwas/tests/test_data_sync.py -v
"""
import json

import pytest

from auric_gwas import data_sync as ds, paths


def test_manifest_exists_and_lists_every_vendored_file():
    assert ds.MANIFEST.is_file(), f"missing {ds.MANIFEST}; run scripts/sync_from_staph_snp_db.py"
    recorded = set(ds.load_manifest().get("files", {}))
    assert recorded == set(ds.VENDORED_FILES), (
        f"manifest files {sorted(recorded)} != expected {sorted(ds.VENDORED_FILES)}")


def test_vendored_files_match_manifest():
    problems = ds.verify_against_manifest()
    assert not problems, "vendored data drifted from PROVENANCE.json:\n  " + "\n  ".join(problems)


def test_source_matches_manifest():
    root = ds.staph_snp_db_root()
    if root is None:
        pytest.skip("staph_snp_db checkout not found ($STAPH_SNP_DB or ../staph_snp_db)")
    problems = ds.verify_against_source(root)
    assert not problems, (
        "live staph_snp_db source drifted from PROVENANCE.json — re-run "
        "scripts/sync_from_staph_snp_db.py:\n  " + "\n  ".join(problems))


# --- runtime bundle checksum assertion -------------------------------------------------------

def _bundle_json(manifest, *, backbone=None, profile_set=None, db_version="9.9.9"):
    """Synthesize a bundle db_version.json; default checksums match the manifest (compatible)."""
    return {
        "db_version": db_version,
        "backbone_panel": {"sha256": backbone or manifest["backbone_panel_sha256"]},
        "coordinate_system": {"profile_hmm_checksums": {
            "profile_set_sha256": profile_set or manifest["coordinate_profile_set_sha256"]}},
    }


def _write(tmp_path, obj):
    p = tmp_path / "db_version.json"
    p.write_text(json.dumps(obj))
    return p


def test_bundle_compatible_passes_on_matching_freeze(tmp_path):
    m = ds.load_manifest()
    paths.assert_bundle_compatible(_write(tmp_path, _bundle_json(m)))  # no raise


def test_bundle_mismatch_on_backbone_raises(tmp_path):
    m = ds.load_manifest()
    p = _write(tmp_path, _bundle_json(m, backbone="deadbeef"))
    with pytest.raises(paths.BundleMismatchError):
        paths.assert_bundle_compatible(p)


def test_bundle_mismatch_on_coordinate_raises(tmp_path):
    m = ds.load_manifest()
    p = _write(tmp_path, _bundle_json(m, profile_set="deadbeef"))
    with pytest.raises(paths.BundleMismatchError):
        paths.assert_bundle_compatible(p)


def test_missing_bundle_raises_filenotfound(tmp_path):
    with pytest.raises(FileNotFoundError):
        paths.assert_bundle_compatible(tmp_path / "absent.json")


def test_skip_env_bypasses(tmp_path, monkeypatch):
    m = ds.load_manifest()
    p = _write(tmp_path, _bundle_json(m, backbone="deadbeef"))
    monkeypatch.setenv("AURIC_SKIP_BUNDLE_CHECK", "1")
    paths.assert_bundle_compatible(p)  # no raise despite mismatch


def test_live_bundle_matches_manifest():
    """If the default/development bundle is mounted, it must match the vendored freeze."""
    bundle = paths.db_version_json()
    if not bundle.is_file():
        pytest.skip(f"no AURIC bundle at {bundle}")
    paths.assert_bundle_compatible(bundle)  # raises BundleMismatchError on drift
