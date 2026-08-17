#!/usr/bin/env python3
"""Re-copy the vendored family-definition files from a staph_snp_db checkout and rewrite the
provenance manifest (auric_gwas/data/PROVENANCE.json).

Run this only after a staph_snp_db panel/coordinate re-freeze, or to seed the manifest the first
time. It resolves the source repo from $STAPH_SNP_DB, else the sibling ../staph_snp_db.

    python scripts/sync_from_staph_snp_db.py           # copy + rewrite manifest
    python scripts/sync_from_staph_snp_db.py --check    # report drift, write nothing (exit 1 if any)

The mapping of vendored file -> source path lives in auric_gwas/data_sync.py, shared with the test.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys

from auric_gwas import data_sync as ds


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="report drift vs the live source and exit non-zero; copy nothing")
    args = ap.parse_args()

    root = ds.staph_snp_db_root()
    if root is None:
        print("error: staph_snp_db not found. Set $STAPH_SNP_DB or place it at ../staph_snp_db.",
              file=sys.stderr)
        return 2
    print(f"staph_snp_db: {root}  (db_version {ds.live_db_version(root)})")

    if args.check:
        if not ds.MANIFEST.is_file():
            print("no manifest yet; run without --check to create it.", file=sys.stderr)
            return 1
        problems = ds.verify_against_manifest() + ds.verify_against_source(root)
        if problems:
            print("DRIFT:", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            return 1
        print("in sync: vendored files match manifest and live source.")
        return 0

    for name, rel in ds.VENDORED_FILES.items():
        src = root / rel
        dst = ds.DATA_DIR / name
        shutil.copyfile(src, dst)
        print(f"  copied {rel} -> data/{name}  ({ds.sha256_file(dst)[:12]}…)")

    manifest = ds.build_manifest(root)
    ds.MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {ds.MANIFEST.relative_to(ds.DATA_DIR.parents[1])}  "
          f"(db_version {manifest['db_version']})")

    problems = ds.verify_against_manifest()
    if problems:
        print("post-write verification FAILED:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
