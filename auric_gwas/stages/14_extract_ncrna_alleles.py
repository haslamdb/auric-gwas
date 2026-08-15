#!/usr/bin/env python3
"""CHOP NICU projection — ncRNA front-end. Fork of scripts/extract_ncrna_alleles_p2.py.

Behaviour is unchanged from the panel script; only the output paths are rebound to
/fastpool/sausnp/chop_nicu_reanalysis/ and the genome list is taken as input instead of
read from results/qc/stage3_input_list.tsv (which is panel state and must not be touched).

What it does, unchanged:
  ncRNAs are dispersed and are not reliably flanked by core genes, so each genome's copy is
  located by scanning the assembly with nhmmer (the DNA-HMM search tool; the match-state
  numbering has to match model.hmm, which is why the profile is an HMM and not a covariance
  model). Per genome: nhmmer a concatenated query of the 46 SINGLE-COPY family HMMs against
  the assembly, keep the best hit per family by E-value, extract and strand-orient the hit
  sequence, append it to buckets_noncoding/<family>.fna. The 10 multi-copy families are
  masked for v1.

Read-only inputs (never written): /fastpool/sausnp/db/coordinate_table.parquet,
/fastpool/sausnp/profiles/ncrna/<fam>/model.hmm, the assemblies.

Differences from the panel script, all deliberate:
  * NC_HMM is built into the COHORT directory, not into /fastpool/sausnp/db/genotype/.
    The frozen panel copy is compared byte-for-byte and the verdict is printed; the frozen
    file is opened read-only and never rewritten.
  * The genome list comes from --genome-list (TSV with genome_id,path) or, by default, from
    globbing the assemblies directory. Assemblies are produced concurrently, so the set is
    snapshotted at start and the snapshot is written to the log directory.
  * Buckets are appended to, as in the panel script, but a per-arm done-list makes a re-run
    resumable without writing a genome into a bucket twice.

Usage:
  14_extract_ncrna_alleles.py [--genome-list T] [--limit N] [--bench N] [--workers 16]
"""
from __future__ import annotations

import argparse
import os
import resource
import statistics
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

# --- paths -------------------------------------------------------------------------------
# Path constants come from scripts/chop_nicu_reanalysis/common_paths.py (owned by another agent).
# The except branch is the identical set, kept so this script still runs standalone.
try:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    from common_paths import (ASSEMBLIES, BUCKETS_NC, COH, COORD, HMMBIN, INPUT, NC_HMM as
                              PANEL_NC_HMM, PREDICT, PROF_NCRNA)
except Exception:
    COH = Path("/fastpool/sausnp/chop_nicu_reanalysis")
    ASSEMBLIES = COH / "assemblies"
    PREDICT = COH / "predict"
    INPUT = COH / "input_list.tsv"
    BUCKETS_NC = COH / "buckets_noncoding"
    COORD = Path("/fastpool/sausnp/db/coordinate_table.parquet")
    PROF_NCRNA = Path("/fastpool/sausnp/profiles/ncrna")
    PANEL_NC_HMM = Path("/fastpool/sausnp/db/genotype/nc_single.hmm")  # read-only
    HMMBIN = "/home/david/miniforge3/envs/abricate/bin"

BUCKETS = BUCKETS_NC
LOGS = COH / "noncoding_logs"          # kept separate from the coding arm's COH/logs
NC_HMM = COH / "nc_single.hmm"         # cohort-owned copy; PANEL_NC_HMM is only ever read
DONE_LIST = LOGS / "ncrna_done.txt"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
EVAL = 1e-5


def revcomp(s):
    return s.translate(COMP)[::-1]


def contigs_by_name(fna):
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


def build_query_hmm():
    """Concatenate the 46 single-copy ncRNA HMMs into the cohort's own query file.

    Identical selection rule to the panel script (feature_id without '__c' == single copy),
    so the file content is expected to be byte-identical to the frozen panel copy. That
    equality is asserted below rather than assumed, and the frozen copy is only read.
    """
    ct = pd.read_parquet(COORD, columns=["feature_id", "feature_type"])
    nc = ct[ct.feature_type == "ncRNA"].feature_id.unique()
    single = sorted(set(f for f in nc if "__c" not in f))
    NC_HMM.parent.mkdir(parents=True, exist_ok=True)
    with open(NC_HMM, "w") as out:
        for fam in single:
            hmm = PROF_NCRNA / fam / "model.hmm"
            if hmm.exists():
                out.write(hmm.read_text())
    if PANEL_NC_HMM.exists():
        same = PANEL_NC_HMM.read_bytes() == NC_HMM.read_bytes()
        print(f"nc_single.hmm identical to frozen panel copy: {same}", flush=True)
    return single


def extract_genome(gid_path):
    """nhmmer the single-copy ncRNA HMMs vs the assembly -> {family: best-hit oriented seq}."""
    gid, gpath = gid_path
    if not Path(gpath).exists():
        return gid, {}
    env = dict(os.environ, PATH=f"{HMMBIN}:{os.environ.get('PATH','')}")
    with tempfile.TemporaryDirectory() as wd:
        tbl = os.path.join(wd, "t.tbl")
        r = subprocess.run([f"{HMMBIN}/nhmmer", "--tblout", tbl, "--cpu", "1", "--noali",
                            "-E", str(EVAL), str(NC_HMM), str(gpath)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        if r.returncode != 0 or not os.path.exists(tbl):
            return gid, {}
        best = {}   # family -> (evalue, target, envfrom, envto, strand)
        for ln in open(tbl):
            if ln.startswith("#"):
                continue
            f = ln.split()
            if len(f) < 14:
                continue
            target, family = f[0], f[2]
            envfrom, envto, strand = int(f[8]), int(f[9]), f[11]
            ev = float(f[12])
            if family not in best or ev < best[family][0]:
                best[family] = (ev, target, envfrom, envto, strand)
    contigs = contigs_by_name(gpath)
    out = {}
    for fam, (ev, target, ef, et, strand) in best.items():
        contig = contigs.get(target)
        if contig is None:
            continue
        lo, hi = (ef, et) if ef <= et else (et, ef)
        seq = contig[lo - 1:hi]
        if strand == "-":
            seq = revcomp(seq)
        if seq:
            out[fam] = seq
    return gid, out


def resolve_genomes(args):
    """(genome_id, assembly_path) pairs, snapshotted now — assemblies are still arriving."""
    src = args.genome_list or (INPUT if Path(INPUT).exists() else None)
    if src:
        df = pd.read_csv(src, sep="\t", dtype=str)
        items = list(zip(df.genome_id.astype(str), df.path.astype(str)))
        items = [it for it in items if Path(it[1]).exists()]
    else:
        items = [(p.stem, str(p)) for p in sorted(ASSEMBLIES.glob("*.fna"))]
    if args.require_predict:
        have = {p.name for p in PREDICT.iterdir()
                if (p / f"{p.name}.assign.tsv").exists()} if PREDICT.exists() else set()
        items = [it for it in items if it[0] in have]
    return items


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-list",
                    help="TSV with genome_id + path; default = COH/input_list.tsv, else glob assemblies/")
    ap.add_argument("--require-predict", action="store_true",
                    help="restrict to genomes that already have predict/<gid>/<gid>.assign.tsv")
    ap.add_argument("--bench", type=int, default=0, help="first N genomes, no bucket writes")
    ap.add_argument("--limit", type=int, default=0, help="first N genomes, buckets written")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(8192, hard), hard))
    LOGS.mkdir(parents=True, exist_ok=True)

    single = build_query_hmm()
    print(f"query: {len(single)} single-copy ncRNA HMMs -> {NC_HMM}", flush=True)

    items = resolve_genomes(args)
    n_snapshot = len(items)
    done = set(DONE_LIST.read_text().split()) if DONE_LIST.exists() else set()
    if done and not args.bench:
        items = [it for it in items if it[0] not in done]
    if args.bench:
        items = items[:args.bench]
    elif args.limit:
        items = items[:args.limit]
    (LOGS / "ncrna_genome_snapshot.txt").write_text("\n".join(g for g, _ in items) + "\n")
    print(f"nhmmer ncRNA extraction for {len(items)} genomes "
          f"({n_snapshot} in snapshot, {len(done)} already done, {args.workers} workers)",
          flush=True)
    if not items:
        return 0

    BUCKETS.mkdir(parents=True, exist_ok=True)
    fh = {}
    per_ok = []
    t0 = time.time()
    with open(DONE_LIST, "a") as dl, ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (gid, out) in enumerate(ex.map(extract_genome, items, chunksize=4), 1):
            per_ok.append(len(out))
            if not args.bench:
                for fam, seq in out.items():
                    h = fh.get(fam)
                    if h is None:
                        h = fh[fam] = open(BUCKETS / f"{fam}.fna", "a")
                    h.write(f">{gid}\n{seq}\n")
                dl.write(gid + "\n")
            if i % 200 == 0:
                print(f"  {i}/{len(items)}  ({(time.time()-t0)/i:.2f} s/genome)", flush=True)
    for h in fh.values():
        h.close()
    el = time.time() - t0
    med = statistics.median(per_ok) if per_ok else 0
    fams = sorted({p.stem for p in BUCKETS.glob("*.fna")} & set(single))
    print(f"\nncRNA: median {med:.0f}/{len(single)} families hit per genome "
          f"({100*med/max(len(single),1):.0f}%)")
    print(f"families with at least one sequence in a bucket: {len(fams)}/{len(single)}")
    print(f"elapsed {el:.1f}s for {len(items)} genomes = {el/max(len(items),1):.2f} s/genome")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
