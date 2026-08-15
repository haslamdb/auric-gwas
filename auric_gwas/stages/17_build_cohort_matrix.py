#!/usr/bin/env python3
"""CHOP NICU projection, step 17 — build the COHORT genotype matrix on the FROZEN panel column order.

This is a fork of scripts/build_genotype_matrix_stage9.py (walk + densify) with the
scripts/v050_narrow/06_matrix_v050.py profiles_masked fallback folded in. Three things differ from
the panel builder, and only three:

  1. ROWS are the cohort genomes (the CHOP NICU SRA accessions that have an assembly and a bucket
     entry at run time), not the 20,685 panel genomes. The count is read from the input list, never
     hardcoded, so the script is re-runnable as the assembler produces more genomes.

  2. COLUMNS are READ from the frozen /fastpool/sausnp/db/genotype_matrix/variant_order.parquet
     (1,231,021 columns) rather than re-derived from variants.parquet. Re-deriving would reproduce
     the same set today, but the whole point of a projection is that cohort column j and panel
     column j are the SAME variant; binding to the frozen file makes that guaranteed rather than
     incidental. variant_order.parquet is copied verbatim into the cohort output directory so the
     cohort matrix is self-describing.

  3. BUCKETS are the cohort buckets under /fastpool/sausnp/chop_nicu_reanalysis/{buckets,
     buckets_noncoding}/ built by 12_genotype_call.py and 15_genotype_noncoding.py. Nothing under
     /fastpool/sausnp/db/ is opened for write; the frozen DB is never opened at all (the coordinate
     table is read from parquet).

Everything else is byte-identical to the panel walk: hmmalign --amino for CDS / --dna for
non-coding, COV_MIN 0.80 divergence gate, premature-stop LoF gate, IDENT_MIN 0.70 switched-IGR
gate, and the same 0 / 1 / -1 encoding — including the sub-feature missingness fix of 2026-08-08,
which both builders carry (see build_genotype_matrix_stage9.py's header for the measurement that
motivated it). A match state with no observed base inside an otherwise-called feature is -1, not 0;
before the fix such a cell was densified to reference, which is how five genomes whose AgrA ORF
stops at codon 226 of 238 came to be recorded as reference at every downstream match state of
group_5710. Coding and non-coding features are contiguous column
blocks of ONE matrix, in the frozen (feature_id, match_state, codon_pos, alt) order — the frozen
order interleaves them alphabetically by feature_id, which is why the densify walks
variant_order's own feature ranges rather than concatenating two matrices.

Masks applied:
  * the 54 permanently masked paralog families are already absent from the frozen column set
    (group_4915 / group_5101 were unmasked at db 0.5.0 and ARE columns — hence the
    profiles_masked/ fallback for their model.hmm);
  * config/excluded_features.txt (11 MSCRAMM adhesins) is checked against the column set and
    reported. Columns are NEVER dropped: dropping one would break the frozen column
    correspondence with the panel matrix. Any excluded-feature column is recorded in
    excluded_feature_columns.parquet so downstream analysis can filter with a boolean mask.

Outputs under /fastpool/sausnp/chop_nicu_reanalysis/genotype_matrix/ (--out-dir to redirect) :
  shards/<feature_id>.npz          per-feature called rows + (row, col) ALT cells + (row, col)
                                   gapped-match-state MISSING cells                [restartable]
  genotypes.zarr                   int8 (n_cohort_genomes x 1,231,021)
  genome_order.parquet             row_index, genome_id (+ CC/ST/manifest metadata)
  variant_order.parquet            verbatim copy of the frozen panel column order
  excluded_feature_columns.parquet col_index / feature_id of any excluded-by-name column
  matrix_stats.json                cell counts, missingness, per-feature coverage

Usage: 17_build_cohort_matrix.py [--phase walk|densify|both] [--workers 24] [--force] [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_paths as P  # noqa: E402

OUTDIR = P.MATRIX
SHARDDIR = OUTDIR / "shards"
BLOCK_COLS = 4096 * 16

# module globals set in build_maps() before the pool forks; inherited copy-on-write by the workers
ROW_OF: dict[str, int] = {}
COL_BY_FEAT: dict[str, dict] = {}     # feature_id -> {(match_state, codon_pos, alt): col_index}
REF_BY_FEAT: dict[str, dict] = {}     # feature_id -> {match_state: ref_allele}   (non-coding only)
FT_OF: dict[str, str] = {}            # feature_id -> CDS / IGR / ncRNA


# ------------------------------------------------------------------------------------------------
# FASTA + hmmalign (identical to build_genotype_matrix_stage9.py)
# ------------------------------------------------------------------------------------------------
def read_fasta(path):
    d, name, seq = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name is not None:
                d[name] = "".join(seq)
            name = ln[1:].split()[0]
            seq = []
        else:
            seq.append(ln.strip())
    if name is not None:
        d[name] = "".join(seq)
    return d


def _hmmalign(hmm, fasta, dna_mode):
    env = P.hmm_env()
    with tempfile.TemporaryDirectory() as wd:
        afa = os.path.join(wd, "a.afa")
        mode = "--dna" if dna_mode else "--amino"
        with open(afa, "w") as af:
            r = subprocess.run([P.HMMALIGN, mode, "--trim", "--outformat", "afa",
                                str(hmm), str(fasta)], stdout=af, stderr=subprocess.DEVNULL, env=env)
        if r.returncode != 0:
            return None
        return read_fasta(afa)


# ------------------------------------------------------------------------------------------------
# Per-feature walk
# ------------------------------------------------------------------------------------------------
def cols_by_match_state(colmap):
    """match_state -> the column indices of every catalogued variant at that match state (all codon
    positions, all ALT alleles). Used to encode a gapped match state as missing across its columns."""
    out: dict[int, list[int]] = {}
    for (m, _cp, _alt), col in colmap.items():
        out.setdefault(m, []).append(col)
    return out


def call_cds(fam):
    """CDS gates: coverage(occupied match states / n_match) < 0.80 -> missing; a stop codon at any
    match state before the last -> LoF -> missing. Otherwise the genome is `called` for the whole
    feature block, its ALT columns are set to 1, and the match states it has no observed base at
    (the up-to-20% the coverage gate tolerates) are set to -1 rather than left at reference."""
    pfa, dfa = P.BUCKETS / f"{fam}.faa", P.BUCKETS / f"{fam}.ffn"
    if not pfa.exists() or not dfa.exists():
        return None, None, None
    aln = _hmmalign(P.cds_profile(fam) / "model.hmm", pfa, dna_mode=False)
    if not aln:
        return None, None, None
    dna = read_fasta(dfa)
    s0 = next(iter(aln.values()))
    matchset = {c for c, ch in enumerate(s0) if ch == "-" or ch.isupper()}
    n_ms = len(matchset)
    if n_ms == 0:
        return None, None, None
    colmap = COL_BY_FEAT.get(fam, {})
    ms_cols = cols_by_match_state(colmap)
    called, arow, acol, mrow, mcol = [], [], [], [], []
    for gid, aseq in aln.items():
        row = ROW_OF.get(gid)
        if row is None:
            continue
        d = dna.get(gid)
        if not d:
            continue
        p = occupied = ms = 0
        ms_to_p = {}
        for c in range(len(aseq)):
            ch = aseq[c]
            if c in matchset:
                ms += 1
                if ch != "-":
                    ms_to_p[ms] = p
                    p += 1
                    occupied += 1
            elif ch != ".":
                p += 1
        if occupied / n_ms < P.COV_MIN:
            continue                                   # divergent -> missing
        lof = False
        alleles = {}
        for m, pp in ms_to_p.items():
            codon = d[3 * pp:3 * pp + 3]
            if len(codon) < 3:
                continue
            if codon.upper() in P.STOPS and m < n_ms:
                lof = True
                break
            for j in range(3):
                alleles[(m, j + 1)] = codon[j].upper()
        if lof:
            continue                                   # LoF -> missing
        called.append(row)
        for (m, cp), base in alleles.items():
            col = colmap.get((m, cp, base))
            if col is not None:
                arow.append(row)
                acol.append(col)
        # sub-feature missingness: a match state with no observed base (an aligner gap, or a codon
        # truncated by the end of the CDS) is NOT reference. Encode every column of it as -1.
        obs = {m for (m, _cp) in alleles}
        for m, cols in ms_cols.items():
            if m not in obs:
                for col in cols:
                    mrow.append(row)
                    mcol.append(col)
    return called, (arow, acol), (mrow, mcol)


def call_noncoding(feat):
    """IGR / ncRNA gates: coverage < 0.80 -> missing; identity-to-reference < 0.70 -> switched IGR
    -> missing (the panel treats a switched IGR as categorical, not as a pile of SNPs). Match states
    with no observed ACGT base inside a called feature are -1, not reference."""
    bucket = P.BUCKETS_NC / f"{feat}.fna"
    if not bucket.exists() or bucket.stat().st_size == 0:
        return None, None, None
    aln = _hmmalign(P.PROFROOT[FT_OF[feat].lower()] / feat / "model.hmm", bucket, dna_mode=True)
    if not aln:
        return None, None, None
    s0 = next(iter(aln.values()))
    match_cols = [c for c, ch in enumerate(s0) if ch == "-" or ch.isupper()]
    n_ms = len(match_cols)
    if n_ms == 0:
        return None, None, None
    colmap = COL_BY_FEAT.get(feat, {})
    ms_cols = cols_by_match_state(colmap)
    ref = REF_BY_FEAT.get(feat, {})
    called, arow, acol, mrow, mcol = [], [], [], [], []
    for gid, aseq in aln.items():
        row = ROW_OF.get(gid)
        if row is None:
            continue
        alleles = {}
        occ = matchref = 0
        for ms, c in enumerate(match_cols, 1):
            ch = aseq[c]
            if ch != "-":
                b = ch.upper()
                if b in P.ACGT:
                    alleles[ms] = b
                    occ += 1
                    if b == ref.get(ms):
                        matchref += 1
        if occ / n_ms < P.COV_MIN:
            continue                                   # divergent -> missing
        if (matchref / occ if occ else 0.0) < P.IDENT_MIN:
            continue                                   # switched -> missing (categorical)
        called.append(row)
        for ms, b in alleles.items():
            col = colmap.get((ms, 0, b))
            if col is not None:
                arow.append(row)
                acol.append(col)
        # sub-feature missingness: a match state with no observed ACGT base is NOT reference.
        for ms, cols in ms_cols.items():
            if ms not in alleles:
                for col in cols:
                    mrow.append(row)
                    mcol.append(col)
    return called, (arow, acol), (mrow, mcol)


def process_feature(fid):
    shard = SHARDDIR / f"{fid}.npz"
    if shard.exists():
        return (fid, "cached", 0)
    called, alt, miss = (call_cds(fid) if FT_OF[fid] == "CDS" else call_noncoding(fid))
    if called is None:
        np.savez_compressed(shard, called=np.empty(0, np.uint16), alt=np.empty((0, 2), np.uint32),
                            miss=np.empty((0, 2), np.uint32))
        return (fid, "empty", 0)
    arow, acol = alt
    mrow, mcol = miss
    np.savez_compressed(shard,
                        called=np.asarray(called, np.uint16),
                        alt=np.asarray(list(zip(arow, acol)), np.uint32).reshape(-1, 2),
                        miss=np.asarray(list(zip(mrow, mcol)), np.uint32).reshape(-1, 2))
    return (fid, len(called), len(arow))


# ------------------------------------------------------------------------------------------------
# Maps
# ------------------------------------------------------------------------------------------------
def cohort_genomes(input_list: Path) -> list[str]:
    with open(input_list, newline="") as fh:
        gids = sorted({r["genome_id"].strip() for r in csv.DictReader(fh, delimiter="\t")
                       if r.get("genome_id", "").strip()})
    if not gids:
        raise SystemExit(f"no genomes in {input_list}")
    assert len(gids) <= 65535, "uint16 row index overflow"
    return gids


def build_maps(input_list: Path):
    """Cohort row order + the FROZEN panel column order + per-feature lookup maps."""
    global ROW_OF, COL_BY_FEAT, REF_BY_FEAT, FT_OF
    genome_ids = cohort_genomes(input_list)
    ROW_OF = {g: i for i, g in enumerate(genome_ids)}

    # THE FROZEN COLUMN ORDER — read, never re-derived.
    p = pd.read_parquet(P.VARIANT_ORDER)
    assert list(p.col_index) == list(range(len(p))), "variant_order col_index is not 0..n-1"
    p["match_state"] = p.match_state.astype(int)
    p["codon_pos"] = p.codon_pos.astype(int)

    for r in p.itertuples(index=False):
        COL_BY_FEAT.setdefault(r.feature_id, {})[(r.match_state, r.codon_pos, r.alt)] = r.col_index

    ct = pd.read_parquet(P.COORD, columns=["feature_id", "feature_type", "match_state",
                                           "codon_pos", "ref_allele"])
    ct = ct[ct.feature_id.isin(set(p.feature_id))]
    for r in ct.itertuples(index=False):
        FT_OF[r.feature_id] = r.feature_type
        if r.feature_type != "CDS":
            # the non-coding identity gate is per match state (codon_pos is always 0 there)
            REF_BY_FEAT.setdefault(r.feature_id, {})[int(r.match_state)] = r.ref_allele
    return genome_ids, p


# ------------------------------------------------------------------------------------------------
# Phases
# ------------------------------------------------------------------------------------------------
def phase_walk(workers, feats):
    SHARDDIR.mkdir(parents=True, exist_ok=True)
    todo = [f for f in feats if not (SHARDDIR / f"{f}.npz").exists()]
    print(f"walk: {len(feats)} features, {len(todo)} to do ({workers} workers)", flush=True)
    t0 = time.time()
    done = tot_called = tot_alt = n_empty = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for fid, nc, na in ex.map(process_feature, todo, chunksize=2):
            done += 1
            if isinstance(nc, int):
                tot_called += nc
                tot_alt += na
            else:
                n_empty += 1
            if done % 400 == 0:
                print(f"  {done}/{len(todo)}  ({time.time()-t0:.0f}s)  cum alt {tot_alt:,}", flush=True)
    print(f"walk done: {len(feats)} shards, {n_empty} empty, {tot_alt:,} alt cells "
          f"({time.time()-t0:.0f}s)", flush=True)


def phase_densify(genome_ids, var_meta):
    """Streamed column-block densify (from 06_matrix_v050.py) — never allocates the dense matrix."""
    import zarr
    from zarr.codecs import BloscCodec
    n_g, n_v = len(genome_ids), len(var_meta)
    frange = var_meta.groupby("feature_id", sort=False)["col_index"].agg(["min", "max"]).sort_values("min")
    print(f"densify (streamed): {n_g} x {n_v:,} int8", flush=True)

    zpath = OUTDIR / "genotypes.zarr"
    if zpath.exists():
        shutil.rmtree(zpath)
    g = zarr.open_group(str(zpath), mode="w")
    comp = BloscCodec(cname="zstd", clevel=5, shuffle="bitshuffle")
    arr = g.create_array("genotypes", shape=(n_g, n_v), chunks=(n_g, 4096), dtype="int8",
                         compressors=[comp], fill_value=-1)

    alt_per_col = np.zeros(n_v, np.int64)
    called_per_col = np.zeros(n_v, np.int64)
    gap_per_col = np.zeros(n_v, np.int64)     # gapped match states inside a CALLED feature
    feat_called = {}
    feat_gap = {}
    buf = np.full((n_g, BLOCK_COLS), -1, np.int8)
    buf_start = buf_end = 0
    t0 = time.time()
    checked_block: dict = {}
    n_blocks_checked = 0

    def flush():
        """Write the buffered column block and re-derive its cell counts from the written cells.

        The bookkeeping control is run on EVERY block, not once: the sub-feature missingness fix
        makes `missing` the sum of two independent terms (genomes uncalled for the whole feature,
        plus gapped match states inside called features), and a gap cell that collided with an ALT
        cell would silently cancel between them. Scanning 1,446 x 65,536 int8 three times costs
        well under a second per block."""
        nonlocal buf_start, buf_end, n_blocks_checked
        if buf_end == buf_start:
            return
        sl = buf[:, : buf_end - buf_start]
        arr[:, buf_start:buf_end] = sl
        sc = {"missing": int((sl == -1).sum()), "alt": int((sl == 1).sum()),
              "ref": int((sl == 0).sum())}
        bk = {"missing": int((n_g - called_per_col[buf_start:buf_end]
                              + gap_per_col[buf_start:buf_end]).sum()),
              "alt": int(alt_per_col[buf_start:buf_end].sum()),
              "ref": int((called_per_col[buf_start:buf_end] - alt_per_col[buf_start:buf_end]
                          - gap_per_col[buf_start:buf_end]).sum())}
        assert sc == bk, (f"densify bookkeeping disagrees with the written cells "
                          f"[cols {buf_start}:{buf_end}]: scanned {sc} vs bookkept {bk}")
        n_blocks_checked += 1
        if not checked_block:
            checked_block.update({"cols": [buf_start, buf_end], "scanned": sc, "bookkept": bk})
        buf[:, : buf_end - buf_start] = -1
        buf_start = buf_end

    for i, fid in enumerate(frange.index, 1):
        c0, c1 = int(frange.loc[fid, "min"]), int(frange.loc[fid, "max"]) + 1
        assert c0 == buf_end, f"frozen column blocks are not contiguous at {fid}: {c0} != {buf_end}"
        if c1 - buf_start > BLOCK_COLS:
            flush()
        shard = SHARDDIR / f"{fid}.npz"
        z = np.load(shard)
        if "miss" not in z.files:
            raise SystemExit(f"shard {fid}.npz predates the sub-feature missingness fix (no `miss` "
                             f"array). Re-run the walk with --force (or point --out-dir at a new "
                             f"directory) so every shard is rebuilt.")
        called, alt, miss = z["called"], z["alt"], z["miss"]
        feat_called[fid] = int(called.size)
        feat_gap[fid] = int(miss.size and miss.shape[0])
        o0, o1 = c0 - buf_start, c1 - buf_start
        if called.size:
            buf[np.asarray(called, np.int64)[:, None], np.arange(o0, o1)] = 0
            called_per_col[c0:c1] = called.size
        if miss.size:
            # gapped match states inside a called feature -> -1, NOT reference. Written before the
            # ALT pass; the two sets are disjoint by construction (a match state either has an
            # observed base or it does not) and flush() re-derives the counts to prove it.
            mc = miss[:, 1].astype(np.int64)
            buf[miss[:, 0].astype(np.int64), mc - buf_start] = -1
            gap_per_col[c0:c1] += np.bincount(mc - c0, minlength=c1 - c0)
        if alt.size:
            ac = alt[:, 1].astype(np.int64)
            buf[alt[:, 0].astype(np.int64), ac - buf_start] = 1
            alt_per_col[c0:c1] += np.bincount(ac - c0, minlength=c1 - c0)
        buf_end = c1
        if i % 500 == 0:
            print(f"  {i}/{len(frange)} features ({time.time()-t0:.0f}s)", flush=True)
    flush()
    assert (called_per_col >= alt_per_col + gap_per_col).all(), \
        "a column has more alt+gap cells than called genomes"
    n_alt = int(alt_per_col.sum())
    n_gap = int(gap_per_col.sum())
    n_ref = int((called_per_col - alt_per_col - gap_per_col).sum())
    n_missing = n_g * n_v - n_alt - n_ref

    g.attrs["encoding"] = "0=called_ref/non-this-alt, 1=called_this_alt, -1=missing"
    g.attrs["n_genomes"] = n_g
    g.attrs["n_variants"] = n_v
    g.attrs["cohort"] = "chop_nicu (She et al. 2026), projected onto SauSNP"
    g.attrs["column_set"] = (f"FROZEN panel order read from {P.VARIANT_ORDER} ({n_v} columns); "
                             "cohort column j is panel column j")
    g.attrs["row_order"] = "genome_order.parquet"
    g.attrs["col_order"] = "variant_order.parquet"
    print(f"  cells: ref(0) {n_ref:,} | alt(1) {n_alt:,} | missing(-1) {n_missing:,} "
          f"({100*n_missing/(n_g*n_v):.3f}% missing, {100*n_alt/(n_g*n_v):.4f}% alt)", flush=True)
    print(f"  missing breakdown: {n_missing-n_gap:,} whole-feature uncalled + {n_gap:,} gapped "
          f"match states inside called features ({100*n_gap/(n_g*n_v):.4f}% of all cells)", flush=True)

    # ---- sidecars -----------------------------------------------------------------------------
    go = pd.DataFrame({"row_index": range(n_g), "genome_id": genome_ids})
    if P.MANIFEST.exists():
        man = pd.read_csv(P.MANIFEST, sep="\t", dtype=str)
        keep = [c for c in ["SRA_Accession", "Genome", "BioSample_Accession", "bioproject",
                            "MLST_ST", "ST", "CC", "MRSA", "patientID", "collection_date",
                            "section", "invasive_status"] if c in man.columns]
        go = go.merge(man[keep].rename(columns={"SRA_Accession": "genome_id",
                                                "ST": "ST_published", "CC": "CC_published"}),
                      on="genome_id", how="left")
    mlst = P.COH / "typing" / "mlst.tsv"
    if mlst.exists():
        m = pd.read_csv(mlst, sep="\t", dtype=str)
        cols = {c: f"{c}_sausnp" for c in ("ST", "CC") if c in m.columns}
        gcol = "genome_id" if "genome_id" in m.columns else m.columns[0]
        go = go.merge(m[[gcol] + list(cols)].rename(columns={gcol: "genome_id", **cols}),
                      on="genome_id", how="left")
    go.to_parquet(OUTDIR / "genome_order.parquet", index=False)
    var_meta.to_parquet(OUTDIR / "variant_order.parquet", index=False)

    excl = P.excluded_feature_names()
    ex_cols = var_meta[var_meta.feature_id.isin(excl)][["col_index", "feature_id"]]
    ex_cols.to_parquet(OUTDIR / "excluded_feature_columns.parquet", index=False)

    per_feat = pd.DataFrame({"feature_id": list(feat_called.keys()),
                             "n_called": list(feat_called.values())})
    per_feat["n_gap_cells"] = per_feat.feature_id.map(feat_gap)
    per_feat["feature_type"] = per_feat.feature_id.map(FT_OF)
    per_feat["n_columns"] = per_feat.feature_id.map(
        var_meta.groupby("feature_id").size().to_dict())
    per_feat.to_parquet(OUTDIR / "feature_coverage.parquet", index=False)

    stats = {
        "n_genomes": n_g, "n_variants": n_v,
        "ref_cells": n_ref, "alt_cells": n_alt, "missing_cells": int(n_missing),
        "missing_whole_feature_cells": int(n_missing - n_gap),
        "missing_gapped_match_state_cells": n_gap,
        "pct_missing": round(100 * n_missing / (n_g * n_v), 4),
        "pct_missing_gapped_match_state": round(100 * n_gap / (n_g * n_v), 4),
        "pct_alt": round(100 * n_alt / (n_g * n_v), 5),
        "n_features": int(len(frange)),
        "n_features_CDS": int(sum(FT_OF.get(f) == "CDS" for f in frange.index)),
        "n_features_noncoding": int(sum(FT_OF.get(f) in ("IGR", "ncRNA") for f in frange.index)),
        "n_features_zero_called": int((per_feat.n_called == 0).sum()),
        "excluded_feature_columns": int(len(ex_cols)),
        "excluded_feature_names_present": sorted(set(ex_cols.feature_id)),
        "densify_block_control": checked_block,
        "densify_blocks_checked": n_blocks_checked,
        "frozen_variant_order": str(P.VARIANT_ORDER),
    }
    (OUTDIR / "matrix_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(f"  sidecars written ({time.time()-t0:.0f}s)", flush=True)
    print(json.dumps({k: v for k, v in stats.items() if k != "densify_block_control"}, indent=2))
    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=["walk", "densify", "both"], default="both")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--input-list", default=str(P.INPUT))
    ap.add_argument("--force", action="store_true",
                    help="delete existing shards first (required when the genome set changes: "
                         "shards store ROW INDICES, which shift when rows are added)")
    ap.add_argument("--out-dir", default=str(P.MATRIX),
                    help="write shards, Zarr and sidecars here instead of the default cohort "
                         "matrix directory; used to build a second matrix without destroying the "
                         "first (must be under the cohort root)")
    args = ap.parse_args()
    _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(8192, hard), hard))

    global OUTDIR, SHARDDIR
    OUTDIR = Path(args.out_dir).resolve()
    assert OUTDIR == P.COH or P.COH in OUTDIR.parents, \
        f"--out-dir must be under the cohort root {P.COH}; refusing to write to {OUTDIR}"
    SHARDDIR = OUTDIR / "shards"

    OUTDIR.mkdir(parents=True, exist_ok=True)
    if args.force and SHARDDIR.exists():
        shutil.rmtree(SHARDDIR)

    genome_ids, var_meta = build_maps(Path(args.input_list))
    feats = sorted(COL_BY_FEAT.keys())
    n_cds = sum(FT_OF.get(f) == "CDS" for f in feats)
    print(f"{len(genome_ids)} cohort genomes x {len(var_meta):,} FROZEN columns over {len(feats)} "
          f"features (CDS {n_cds}, noncoding {len(feats)-n_cds}) -> {OUTDIR}", flush=True)
    missing_ft = [f for f in feats if f not in FT_OF]
    assert not missing_ft, f"features with no coordinate feature_type: {missing_ft[:5]}"
    for f in P.UNMASKED_AT_050:
        assert f in COL_BY_FEAT, f"{f} is not in the frozen column set"
        assert (P.cds_profile(f) / "model.hmm").exists(), f"no model.hmm resolvable for {f}"

    # A shard's row indices are only valid for the genome set that produced them.
    stamp = OUTDIR / "shards_genome_set.txt"
    cur = "\n".join(genome_ids)
    if SHARDDIR.exists() and stamp.exists() and stamp.read_text() != cur:
        raise SystemExit("existing shards were built for a DIFFERENT genome set; re-run with --force")

    if args.phase in ("walk", "both"):
        phase_walk(args.workers, feats)
        stamp.write_text(cur)
    if args.phase in ("densify", "both"):
        phase_densify(genome_ids, var_meta)
    print("cohort genotype matrix complete. Nothing under /fastpool/sausnp/db was written.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
