#!/usr/bin/env python3
"""Portable path/binary constants for the vendored AURIC genotyping stages (11-17).

Derived from the validated AURIC cohort-genotyping common_paths, made distributable:
  * every frozen READ-ONLY input resolves under AURIC_HOME (env var, default /fastpool/sausnp) with
    the layout db/ (sausnp.duckdb, coordinate_table.parquet, genotype_matrix/, hmm/, parquet/) and
    profiles/ + profiles_masked/ — the same layout auric_gwas.paths uses;
  * tool binaries resolve from PATH via shutil.which (override with AURIC_PRODIGAL / AURIC_HMMSEARCH /
    AURIC_HMMALIGN / AURIC_HMMPRESS / AURIC_NHMMER);
  * the two small family-definition files (stage4_core_cds.tsv, excluded_features.txt) are vendored in
    auric_gwas/data/.

WRITE targets (COH and everything under it) are PLACEHOLDERS here: auric_gwas.genotype generates a
per-request module that inherits this one and rebinds them to the request root. Nothing under AURIC_HOME
is ever opened for write. Standard library only, no side effects on import (call ensure_dirs()).
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

# --------------------------------------------------------------------------------------------
# Reference location (READ-ONLY). AURIC_HOME points at the unpacked reference bundle.
# --------------------------------------------------------------------------------------------
AURIC_HOME = Path(os.environ.get("AURIC_HOME", "/fastpool/sausnp"))
DBDIR = AURIC_HOME / "db"
DUCKDB = DBDIR / "sausnp.duckdb"          # open ONLY with duckdb.connect(str(DUCKDB), read_only=True)
COORD = DBDIR / "coordinate_table.parquet"
VARIANTS_PQ = DBDIR / "parquet" / "variants.parquet"
FEATURES_PQ = DBDIR / "parquet" / "features.parquet"
CALLABLE_MASK_PQ = DBDIR / "parquet" / "callable_mask.parquet"

PANEL_GENOTYPE = DBDIR / "genotype"
PANEL_MATRIX = DBDIR / "genotype_matrix"
VARIANT_ORDER = PANEL_MATRIX / "variant_order.parquet"   # frozen 1,231,021-column order
GENOME_ORDER = PANEL_MATRIX / "genome_order.parquet"
PANEL_ZARR = PANEL_MATRIX / "genotypes.zarr"

PROF_CDS = AURIC_HOME / "profiles" / "cds"
PROF_CDS_MASKED = AURIC_HOME / "profiles_masked" / "cds"
PROF_IGR = AURIC_HOME / "profiles" / "igr"
PROF_NCRNA = AURIC_HOME / "profiles" / "ncrna"
PROFROOT = {"igr": PROF_IGR, "ncrna": PROF_NCRNA}

CDS_DB = DBDIR / "hmm" / "cds_all.hmm"     # 1,749 models (group_4915/group_5101 supplied separately)
IGR_DB = DBDIR / "hmm" / "igr_all.hmm"
NCRNA_DB = DBDIR / "hmm" / "ncrna_all.hmm"
NC_HMM = PANEL_GENOTYPE / "nc_single.hmm"  # read-only reference; stage 14 builds its own under COH

UNMASKED_AT_050 = ("group_4915", "group_5101")

# Vendored family-definition files (part of the tool, not the bundle).
_DATA = Path(__file__).resolve().parent.parent / "data"
CORE_CDS_TSV = _DATA / "stage4_core_cds.tsv"          # is_core_cds == True -> 1,751 families
EXCLUDED_FEATURES = _DATA / "excluded_features.txt"   # 11 MSCRAMM adhesins, by name
PANEL_QC = _DATA                                      # stage 13 reads stage4_core_igr.tsv from here

# --------------------------------------------------------------------------------------------
# WRITE targets — PLACEHOLDERS. auric_gwas.genotype rebinds all of these per request.
# --------------------------------------------------------------------------------------------
COH = Path(os.environ.get("AURIC_COH", "/tmp/auric_gwas_cohort"))
ASSEMBLIES = COH / "assemblies"
INPUT = COH / "input_list.tsv"
PREDICT = COH / "predict"
BUCKETS = COH / "buckets"
CALLS = COH / "calls"
BUCKETS_NC = COH / "buckets_noncoding"
CALLS_NC = COH / "calls_noncoding"
SWITCHED_NC = COH / "switched_noncoding"
MATRIX = COH / "genotype_matrix"
QC = COH / "qc"
LOGS = COH / "logs"
CDS_DB_UNMASKED = COH / "hmm" / "cds_unmasked050.hmm"
MANIFEST = COH / "__no_manifest__.tsv"    # absent -> stage 17 skips cohort-metadata merge
RESULTS = COH
COHORT_LABEL = "user_cohort"
SOURCE_LABEL = "user"

# --------------------------------------------------------------------------------------------
# Binaries — from PATH (environment.yml provides them), overridable by env.
# --------------------------------------------------------------------------------------------
def _bin(name: str, env: str) -> str:
    return os.environ.get(env) or shutil.which(name) or name  # bare name fails loudly if absent

PRODIGAL = _bin("prodigal", "AURIC_PRODIGAL")
HMMSEARCH = _bin("hmmsearch", "AURIC_HMMSEARCH")
HMMALIGN = _bin("hmmalign", "AURIC_HMMALIGN")
HMMPRESS = _bin("hmmpress", "AURIC_HMMPRESS")
NHMMER = _bin("nhmmer", "AURIC_NHMMER")
HMMBIN = str(Path(HMMALIGN).parent) if Path(HMMALIGN).parent != Path(".") else ""

# --------------------------------------------------------------------------------------------
# Gates and encodings — identical to the panel pipeline; never retune for a cohort.
# --------------------------------------------------------------------------------------------
EVALUE = "1e-5"
COV_MIN = 0.80
IDENT_MIN = 0.70
STOPS = frozenset({"TAA", "TAG", "TGA"})
ACGT = frozenset("ACGT")
GT_REF, GT_ALT, GT_MISSING = 0, 1, -1


def hmm_env() -> dict:
    """os.environ with the HMMER bin directory prepended to PATH (hmmalign calls esl-* helpers)."""
    if HMMBIN:
        return dict(os.environ, PATH=f"{HMMBIN}:{os.environ.get('PATH', '')}")
    return dict(os.environ)


def cds_profile(fam: str) -> Path:
    """Directory holding a CDS family's model.hmm, with the profiles_masked/ fallback db 0.5.0 needs."""
    p = PROF_CDS / fam
    return p if (p / "model.hmm").exists() else PROF_CDS_MASKED / fam


class ProfileDirWithMaskedFallback:
    def __truediv__(self, fam: str) -> Path:
        return cds_profile(fam)

    def __fspath__(self) -> str:
        return str(PROF_CDS)

    def __repr__(self) -> str:
        return f"ProfileDirWithMaskedFallback({PROF_CDS})"


PROF_CDS_FALLBACK = ProfileDirWithMaskedFallback()


def core_cds_families() -> list[str]:
    """The 1,751 genotypable CDS families in stage4_core_cds.tsv order that resolve to a model.hmm."""
    fams, hdr = [], None
    with open(CORE_CDS_TSV) as fh:
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if hdr is None:
                hdr = c
                continue
            row = dict(zip(hdr, c))
            if row.get("is_core_cds", "").strip().lower() in ("true", "1"):
                fams.append(row["gene"])
    return [f for f in fams if (cds_profile(f) / "model.hmm").exists()]


def excluded_feature_names() -> set[str]:
    out = set()
    for ln in open(EXCLUDED_FEATURES):
        ln = ln.split("#", 1)[0].strip()
        if ln:
            out.add(ln)
    return out


def ensure_dirs() -> None:
    for d in (PREDICT, BUCKETS, CALLS, BUCKETS_NC, CALLS_NC, SWITCHED_NC, MATRIX, QC, LOGS,
              CDS_DB_UNMASKED.parent):
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    print(f"AURIC_HOME {AURIC_HOME}")
    print(f"DUCKDB     {DUCKDB}  exists={DUCKDB.exists()}")
    print(f"PROF_CDS   {PROF_CDS}  exists={PROF_CDS.exists()}")
    print(f"CDS_DB     {CDS_DB}  exists={CDS_DB.exists()}")
    print(f"PRODIGAL   {PRODIGAL}")
    print(f"HMMSEARCH  {HMMSEARCH}")
    print(f"core CDS families genotypable: {len(core_cds_families())}")
