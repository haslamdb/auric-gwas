#!/usr/bin/env python3
"""CHOP NICU projection, step 12 — fork of scripts/genotype_call.py with cohort-rebound paths.

Stage 6 phase 3, unchanged algorithm: sweep the cohort predict output into per-family protein and
coding-DNA FASTAs, hmmalign each family's copies to its profile HMM, walk the match states to
codons, and tally per-site nucleotide alleles against the frozen reference allele from
/fastpool/sausnp/db/coordinate_table.parquet (read-only).

Gates preserved exactly as in scripts/genotype_call.py:
  hmmalign --amino --trim
  match-state coverage < 0.80  -> divergent_allele, genome dropped for that family
  a stop codon at any match state before the last -> LoF, genome dropped for that family

Differences from the panel script, all path/scope only:
 1. PREDICT / BUCKETS / CALLS are the cohort directories from common_paths; nothing is written
    under /fastpool/sausnp/db/. The coordinate table and the profile HMMs are read-only inputs.
 2. The family list comes from common_paths.core_cds_families(), which resolves each family's
    model.hmm through the profiles_masked/ fallback. The panel script tests
    (PROFILES / f / 'model.hmm').exists() against profiles/cds only, which silently drops
    group_4915 and group_5101 — the two families db 0.5.0 unmasked. Result: 1,751 families here
    vs the 1,749 that a literal re-run of the panel script would produce today.
 3. bucket() enumerates genomes from the cohort input list rather than from every subdirectory of
    PREDICT, so a half-written predict directory cannot enter a bucket.

Outputs, per family, under <COH>/calls/<fam>.tsv:
  match_state, codon_pos, ref_allele, allele_counts(JSON), n_called, n_lof, n_divergent
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_paths as P  # noqa: E402

STOPS = set(P.STOPS)
COV_MIN = P.COV_MIN


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


def genome_ids(input_list: Path) -> list[str]:
    ids = []
    with open(input_list) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi = hdr.index("genome_id")
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) > gi:
                ids.append(c[gi])
    return ids


def bucket(input_list: Path):
    """Sweep cohort genomes -> per-family prot/dna FASTAs (genome_id as record name)."""
    P.BUCKETS.mkdir(parents=True, exist_ok=True)
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(8192, hard), hard))
    fams = P.core_cds_families()
    ph = {f: open(P.BUCKETS / f"{f}.faa", "w") for f in fams}
    dh = {f: open(P.BUCKETS / f"{f}.ffn", "w") for f in fams}
    famset = set(fams)
    genomes = [g for g in genome_ids(input_list)
               if (P.PREDICT / g / f"{g}.assign.tsv").exists()]
    print(f"bucketing {len(genomes)} genomes into {len(fams)} families", flush=True)
    for i, gid in enumerate(genomes, 1):
        base = P.PREDICT / gid
        prot = read_fasta(base / f"{gid}.faa")
        dna = read_fasta(base / f"{gid}.ffn")
        for ln in open(base / f"{gid}.assign.tsv"):
            if ln.startswith("family"):
                continue
            fam, pid = ln.split("\t")[:2]
            if fam in famset and pid in prot and pid in dna:
                ph[fam].write(f">{gid}\n{prot[pid]}\n")
                dh[fam].write(f">{gid}\n{dna[pid]}\n")
        if i % 200 == 0:
            print(f"  bucketed {i}/{len(genomes)}", flush=True)
    for f in ph.values():
        f.close()
    for f in dh.values():
        f.close()
    return fams, len(genomes)


def match_columns(seqs):
    """Match columns of an hmmalign afa. HMMER makes each column consistently match
    (UPPER/'-') or insert (lower/'.') across all rows, so one row classifies them all."""
    s0 = next(iter(seqs.values()))
    return [c for c, ch in enumerate(s0) if ch == "-" or ch.isupper()]


def call_family(fam):
    out = P.CALLS / f"{fam}.tsv"
    if out.exists():
        return (fam, "cached")
    pfa, dfa = P.BUCKETS / f"{fam}.faa", P.BUCKETS / f"{fam}.ffn"
    if not pfa.exists() or pfa.stat().st_size == 0:
        return (fam, "empty")
    env = P.hmm_env()
    with tempfile.TemporaryDirectory() as wd:
        afa = os.path.join(wd, "a.afa")
        with open(afa, "w") as af:
            r = subprocess.run([P.HMMALIGN, "--amino", "--trim", "--outformat", "afa",
                                str(P.cds_profile(fam) / "model.hmm"), str(pfa)],
                               stdout=af, stderr=subprocess.DEVNULL, env=env)
        if r.returncode != 0:
            return (fam, "hmmalign_fail")
        aln = read_fasta(afa)
    dna = read_fasta(dfa)
    if not aln:
        return (fam, "empty_aln")
    mcols = match_columns(aln)
    n_ms = len(mcols)
    ref = {}
    ct = pd.read_parquet(P.COORD, filters=[("feature_id", "==", fam)])
    for r in ct.itertuples(index=False):
        ref[(r.match_state, r.codon_pos)] = r.ref_allele

    site_counts = defaultdict(Counter)
    n_lof = n_div = n_called = 0
    for gid, aseq in aln.items():
        d = dna.get(gid)
        if not d:
            continue
        p, occupied, ms_to_p = 0, 0, {}
        ms = 0
        matchset = set(mcols)
        for c in range(len(aseq)):
            ch = aseq[c]
            if c in matchset:
                ms += 1
                if ch != "-":
                    ms_to_p[ms] = p; p += 1; occupied += 1
            else:
                if ch != ".":
                    p += 1
        if n_ms == 0 or occupied / n_ms < COV_MIN:
            n_div += 1
            continue
        lof = False
        alleles = {}
        for m, pp in ms_to_p.items():
            codon = d[3 * pp:3 * pp + 3]
            if len(codon) < 3:
                continue
            if codon.upper() in STOPS and m < n_ms:
                lof = True; break
            for j in range(3):
                alleles[(m, j + 1)] = codon[j].upper()
        if lof:
            n_lof += 1
            continue
        n_called += 1
        for site, base in alleles.items():
            site_counts[site][base] += 1

    with open(out, "w") as fh:
        fh.write("match_state\tcodon_pos\tref_allele\tallele_counts\tn_called\tn_lof\tn_divergent\n")
        for (m, cp), counts in sorted(site_counts.items()):
            ra = ref.get((m, cp), "")
            fh.write(f"{m}\t{cp}\t{ra}\t{json.dumps(dict(counts))}\t{n_called}\t{n_lof}\t{n_div}\n")
    return (fam, "ok")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(P.INPUT))
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--skip-bucket", action="store_true")
    args = ap.parse_args()
    P.CALLS.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    if args.skip_bucket:
        fams, ngen = [p.stem for p in P.BUCKETS.glob("*.faa")], -1
    else:
        fams, ngen = bucket(Path(args.input))
    tb = time.perf_counter() - t0
    print(f"bucketing took {tb:.1f}s", flush=True)
    print(f"calling {len(fams)} families ({args.workers} workers)", flush=True)
    t1 = time.perf_counter()
    ok = 0
    bad = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (fam, st) in enumerate(ex.map(call_family, fams), 1):
            ok += st in ("ok", "cached")
            if st not in ("ok", "cached"):
                bad.append((fam, st))
                print(f"  {fam}: {st}", flush=True)
            if i % 200 == 0:
                print(f"  called {i}/{len(fams)}", flush=True)
    tc = time.perf_counter() - t1
    print(f"done: {ok}/{len(fams)} families -> {P.CALLS}")
    print(f"calling took {tc:.1f}s; total {tb + tc:.1f}s for {ngen} genomes")
    if bad:
        print(f"non-ok families: {bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
