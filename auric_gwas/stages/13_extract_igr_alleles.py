#!/usr/bin/env python3
"""CHOP NICU projection — IGR front-end. Fork of scripts/extract_igr_alleles_p2.py.

Behaviour is unchanged from the panel script; only the output paths are rebound to
/fastpool/sausnp/chop_nicu_reanalysis/ and the genome list is taken as input instead of
read from results/qc/stage3_input_list.tsv (panel state, must not be touched).

What it does, unchanged. This is anchor-and-align, not a Piggy re-run. Piggy re-clustering
over the cohort would produce IGR clusters incompatible with the frozen 905-cluster
coordinate system, so instead each genome's copy of a cluster is sliced out of its assembly
between the two CORE CDS families that flank the cluster in the reference:

  1. flank index: cluster -> (upstream_family, downstream_family), read from the NCTC 8325
     (USA300 fallback) member header in profiles/igr/<C>/aln.fasta and mapped locus-tag ->
     feature_id through the features table. Only clusters whose BOTH flanks are core families
     are anchorable; the rest emit no genotype, which reduces AN exactly as an absent CDS does.
  2. per genome: from predict/<gid>/{.ffn header coordinates, .assign.tsv family->protein_id}
     locate the two flank genes. Require same contig, adjacency, and gap < GAP_MAX = 1500 bp.
     Slice the intergenic region from the assembly and ORIENT it against the reference IGR by
     12-mer overlap, because hmmalign --dna does not search the reverse strand. Append to
     buckets_noncoding/<cluster>.fna as '>genome_id'.

Read-only inputs (never written): /fastpool/sausnp/db/coordinate_table.parquet,
/fastpool/sausnp/db/parquet/features.parquet, results/qc/stage4_core_igr.tsv,
/fastpool/sausnp/profiles/igr/, the assemblies, and predict/<gid>/.

Differences from the panel script, all deliberate:
  * PREDICT, BUCKETS point into the cohort directory.
  * The genome list comes from --genome-list (TSV genome_id,path) or, by default, from the
    intersection of predict/ and the assemblies directory, snapshotted at start because
    assemblies are still being produced.
  * --predict-dir overrides PREDICT so the code path can be validated against a scratch
    predict tree without colliding with the coding arm's output.
  * A per-arm done-list makes a re-run resumable without appending a genome twice.
  * The open-file limit is raised, because ~700 bucket handles are held open at once.

Usage:
  13_extract_igr_alleles.py [--genome-list T] [--predict-dir D] [--bench N] [--limit N]
                            [--workers 16]
"""
from __future__ import annotations

import argparse
import re
import resource
import statistics
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

# --- paths -------------------------------------------------------------------------------
# Path constants come from scripts/chop_nicu_reanalysis/common_paths.py (owned by another agent).
# The except branch is the identical set, kept so this script still runs standalone.
try:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from common_paths import (ASSEMBLIES, BUCKETS_NC, COH, COORD, FEATURES_PQ, INPUT,
                              PANEL_QC, PREDICT, PROF_IGR)
except Exception:
    COH = Path("/fastpool/sausnp/chop_nicu_reanalysis")
    ASSEMBLIES = COH / "assemblies"
    PREDICT = COH / "predict"
    INPUT = COH / "input_list.tsv"
    BUCKETS_NC = COH / "buckets_noncoding"
    COORD = Path("/fastpool/sausnp/db/coordinate_table.parquet")
    FEATURES_PQ = Path("/fastpool/sausnp/db/parquet/features.parquet")
    PROF_IGR = Path("/fastpool/sausnp/profiles/igr")
    PANEL_QC = Path("/home/david/projects/staph_snp_db/results/qc")

FEATURES = FEATURES_PQ                 # read-only
CORE_IGR_TSV = PANEL_QC / "stage4_core_igr.tsv"   # read-only panel state
BUCKETS = BUCKETS_NC
LOGS = COH / "noncoding_logs"          # kept separate from the coding arm's COH/logs
DONE_LIST = LOGS / "igr_done.txt"

GAP_MAX = 1500
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
K = 12

FLANK: dict = {}       # cluster -> (fam_up, fam_down)  (feature_ids; order = reference 5'->3')
REFSEQ: dict = {}      # cluster -> reference IGR sequence (canonical orientation)
REFKM: dict = {}       # cluster -> k-mer set of REFSEQ (precomputed once, for orientation)


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    d, name, seq = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if name is not None:
                d[name] = "".join(seq)
            name = ln[1:].split()[0]; seq = []
        else:
            seq.append(ln.strip())
    if name is not None:
        d[name] = "".join(seq)
    return d


def contigs_by_name(fna):
    """{contig_name (fasta header first token): sequence}. Prodigal seqname == this token, so the
    per-gene ID '<contig_name>_<gene_ordinal>' resolves by rsplit('_',1)[0] for BOTH integer-named
    (in-house '1') and accession-named (NCBI 'NC_007622.1') contigs, and for the SPAdes contigs
    of this cohort, which are renamed '<SRR>_NODE_n_length_...' by 02_assemble.sh."""
    out, name, seq = {}, None, []
    for ln in open(fna):
        if ln.startswith(">"):
            if name is not None:
                out[name] = "".join(seq).upper()
            name = ln[1:].split()[0]; seq = []
        else:
            seq.append(ln.strip())
    if name is not None:
        out[name] = "".join(seq).upper()
    return out


def build_flank_index():
    """cluster -> (fam_up, fam_down) using the reference member's flank locus tags + features map."""
    fe = pd.read_parquet(FEATURES, columns=["feature_id", "nctc8325_locus_tag", "usa300_locus_tag",
                                            "nctc8325_start"])
    nctc_lt2f = {lt: f for lt, f in zip(fe.nctc8325_locus_tag, fe.feature_id) if lt}
    usa_lt2f = {lt: f for lt, f in zip(fe.usa300_locus_tag, fe.feature_id) if lt}
    start_of = dict(zip(fe.feature_id, fe.nctc8325_start))
    igr = pd.read_csv(CORE_IGR_TSV, sep="\t")
    idcol = igr.columns[0]
    clusters = [c for c in igr.loc[igr.is_core_igr == True, idcol].astype(str)
                if (PROF_IGR / c / "model.hmm").exists()]
    flank = {}
    for c in clusters:
        aln = PROF_IGR / c / "aln.fasta"
        fam_pair = None
        for gid_pref, lt2f in (("GCF_000013425.1", nctc_lt2f), ("GCF_000013465.1", usa_lt2f)):
            hdr = next((h for h in read_fasta(aln) if h.startswith(gid_pref + "_")), None)
            if not hdr:
                continue
            tags = re.findall(r"[A-Z]+_\d+", hdr.split(gid_pref + "_", 1)[1])
            fams = [lt2f.get(t) for t in tags[:2]]
            if len(fams) == 2 and all(fams):
                fam_pair = fams
                break
        if not fam_pair:
            continue
        # order by reference start (5'->3'); if unknown keep as-is
        fa, fb = fam_pair
        sa, sb = start_of.get(fa), start_of.get(fb)
        if sa is not None and sb is not None and sa > sb:
            fa, fb = fb, fa
        flank[c] = (fa, fb)
    # reference IGR sequence per cluster (canonical orientation), from coordinate_table ref_allele
    ct = pd.read_parquet(COORD, columns=["feature_id", "feature_type", "match_state", "ref_allele"])
    ct = ct[ct.feature_type == "IGR"].sort_values(["feature_id", "match_state"])
    refseq = {f: "".join(g.ref_allele) for f, g in ct.groupby("feature_id")}
    refkm = {f: kmers(s) for f, s in refseq.items()}
    return flank, refseq, refkm


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)} if len(s) >= K else set()


def orient(seq, refkm):
    if not refkm or len(seq) < K:
        return seq
    km = kmers(seq)
    f = len(km & refkm)
    r = len(kmers(revcomp(seq)) & refkm)
    return seq if f >= r else revcomp(seq)


def extract_genome(gid_path):
    """Return (gid, {cluster: oriented_igr_seq}, status_counter)."""
    gid, gpath = gid_path
    base = PREDICT / gid
    st = Counter()
    if not (base / f"{gid}.assign.tsv").exists() or not Path(gpath).exists():
        return gid, {}, Counter({"no_input": 1})
    # family -> protein_id
    fam2pid = {}
    for ln in open(base / f"{gid}.assign.tsv"):
        if ln.startswith("family"):
            continue
        f = ln.split("\t")
        fam2pid[f[0]] = f[1]
    # protein_id -> (contig_idx, start, end)  from .ffn headers ('>C_g # start # end # strand # ...')
    pid_coord = {}
    for ln in open(base / f"{gid}.ffn"):
        if not ln.startswith(">"):
            continue
        parts = ln[1:].split("#")
        pid = parts[0].strip()
        try:
            start, end = int(parts[1]), int(parts[2])
        except (IndexError, ValueError):
            continue
        cname = pid.rsplit("_", 1)[0]                  # contig name (matches the .fna header token)
        pid_coord[pid] = (cname, start, end)
    contigs = contigs_by_name(gpath)

    out = {}
    for c, (fa, fb) in FLANK.items():
        pa, pb = fam2pid.get(fa), fam2pid.get(fb)
        if not pa or not pb:
            st["missing_family"] += 1; continue
        ca, cb = pid_coord.get(pa), pid_coord.get(pb)
        if not ca or not cb:
            st["no_coord"] += 1; continue
        if ca[0] != cb[0]:
            st["diff_contig"] += 1; continue
        contig = contigs.get(ca[0])
        if contig is None:
            st["bad_contig"] += 1; continue
        lo_end, hi_start = (ca[2], cb[1]) if ca[1] < cb[1] else (cb[2], ca[1])
        s0, s1 = lo_end, hi_start - 1                 # 1-based inclusive between-gene region
        if s1 - s0 + 1 <= 0 or s1 - s0 + 1 > GAP_MAX:
            st["gap"] += 1; continue
        seq = contig[s0:s1 + 1]                         # 0-based slice [lo_end .. hi_start-2]
        if not seq or set(seq) - set("ACGTN"):
            seq = re.sub(r"[^ACGTN]", "N", seq)
        seq = orient(seq, REFKM.get(c))
        out[c] = seq
        st["ok"] += 1
    return gid, out, st


def resolve_genomes(args):
    """(genome_id, assembly_path) pairs, snapshotted now — assemblies are still arriving."""
    src = args.genome_list or (INPUT if Path(INPUT).exists() else None)
    if src:
        df = pd.read_csv(src, sep="\t", dtype=str)
        items = list(zip(df.genome_id.astype(str), df.path.astype(str)))
        items = [it for it in items if Path(it[1]).exists()]
    else:
        items = [(p.stem, str(p)) for p in sorted(ASSEMBLIES.glob("*.fna"))]
    have = {p.name for p in PREDICT.iterdir()
            if (p / f"{p.name}.assign.tsv").exists()} if PREDICT.exists() else set()
    return [it for it in items if it[0] in have]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-list",
                    help="TSV with genome_id + path; default = COH/input_list.tsv, else glob assemblies/")
    ap.add_argument("--predict-dir", help="override the predict/ tree (validation runs)")
    ap.add_argument("--buckets-dir", help="override buckets_noncoding/ (validation runs)")
    ap.add_argument("--bench", type=int, default=0, help="first N genomes, no bucket writes")
    ap.add_argument("--limit", type=int, default=0, help="first N genomes, buckets written")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    global FLANK, REFSEQ, REFKM, PREDICT, BUCKETS, DONE_LIST
    if args.predict_dir:
        PREDICT = Path(args.predict_dir)
    if args.buckets_dir:
        BUCKETS = Path(args.buckets_dir)
        DONE_LIST = BUCKETS.parent / "igr_done.txt"
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(8192, hard), hard))
    LOGS.mkdir(parents=True, exist_ok=True)

    t_idx = time.time()
    FLANK, REFSEQ, REFKM = build_flank_index()
    print(f"flank index: {len(FLANK)}/905 IGR clusters anchorable (both flanks core); "
          f"{len(REFSEQ)} ref sequences  [{time.time()-t_idx:.1f}s]", flush=True)

    items = resolve_genomes(args)
    n_snapshot = len(items)
    done = set(DONE_LIST.read_text().split()) if DONE_LIST.exists() else set()
    if done and not args.bench:
        items = [it for it in items if it[0] not in done]
    if args.bench:
        items = items[:args.bench]
    elif args.limit:
        items = items[:args.limit]
    (LOGS / "igr_genome_snapshot.txt").write_text("\n".join(g for g, _ in items) + "\n")
    print(f"extracting IGRs for {len(items)} genomes "
          f"({n_snapshot} with predict output, {len(done)} already done, {args.workers} workers)",
          flush=True)
    if not items:
        return 0

    BUCKETS.mkdir(parents=True, exist_ok=True)
    bucket_fh = {}
    total = Counter()
    per_genome_ok = []
    t0 = time.time()
    with open(DONE_LIST, "a") as dl, ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (gid, out, st) in enumerate(ex.map(extract_genome, items, chunksize=8), 1):
            total.update(st)
            per_genome_ok.append(st.get("ok", 0))
            if not args.bench:                          # only write buckets in a real run
                for c, seq in out.items():
                    fh = bucket_fh.get(c)
                    if fh is None:
                        fh = bucket_fh[c] = open(BUCKETS / f"{c}.fna", "a")
                    fh.write(f">{gid}\n{seq}\n")
                dl.write(gid + "\n")
            if i % 200 == 0:
                print(f"  {i}/{len(items)} genomes  ({(time.time()-t0)/i:.2f} s/genome)", flush=True)
    for fh in bucket_fh.values():
        fh.close()
    el = time.time() - t0
    med = statistics.median(per_genome_ok) if per_genome_ok else 0
    print(f"\nextracted: median {med:.0f}/{len(FLANK)} IGRs per genome "
          f"({100*med/max(len(FLANK),1):.0f}% of anchorable)")
    print(f"status totals: {dict(total)}")
    print(f"elapsed {el:.1f}s for {len(items)} genomes = {el/max(len(items),1):.2f} s/genome")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
