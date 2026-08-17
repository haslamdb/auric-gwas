"""Guard the three vendored family-definition files against silent drift from staph_snp_db.

The manifest check (test_vendored_files_match_manifest) is self-contained and always runs.
The source check (test_source_matches_manifest) skips unless a staph_snp_db checkout is present.

Run: python -m pytest auric_gwas/tests/test_data_sync.py -v
"""
import pytest

from auric_gwas import data_sync as ds


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
