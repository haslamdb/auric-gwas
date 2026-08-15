#!/usr/bin/env python3
"""CHOP NICU projection, step 11 — fork of scripts/genotype_predict.py with cohort-rebound paths.

Identical algorithm and identical gates to the panel's Stage 6 phase 1-2:
  prodigal -p single  ->  <gid>.faa / <gid>.ffn
  hmmsearch -E 1e-5 --cpu 2 against the CDS profile database
  per family, keep the single best-scoring protein hit; hmm_cov = (hmm_to - hmm_from + 1) / hmm_len

Two deliberate differences from scripts/genotype_predict.py, both required to reproduce the db
0.5.0 panel state rather than the db 0.1.0 state the original script was written against:

 1. Path rebinding. Every write target comes from common_paths (PREDICT under
    /fastpool/sausnp/chop_nicu_reanalysis/). Nothing is written under /fastpool/sausnp/db/.
    The genome list comes from --input (default common_paths.INPUT), never from
    results/qc/stage3_input_list.tsv.

 2. Second scan for the two families db 0.5.0 unmasked. /fastpool/sausnp/db/hmm/cds_all.hmm holds
    1,749 pressed models and does NOT contain group_4915 or group_5101, whose model.hmm files
    remain under profiles_masked/cds/. The panel's own assign.tsv files predate the 0.2.0 mask and
    do carry those two families (1,805 rows). To land in the same place, this script runs a second
    hmmsearch of the same .faa against a two-model HMM file and merges the best hits into the same
    assign.tsv. E-values are unaffected by splitting the scan: hmmsearch scales the E-value by the
    size of the *target sequence* database (the genome's .faa), which is byte-identical in both
    passes, and the "best hit" reduction is per-family so it cannot mix across the two passes.

Outputs, per genome, under <COH>/predict/<gid>/:
  <gid>.faa, <gid>.ffn, <gid>.assign.tsv (family, protein_id, hmm_cov, score), <gid>.unmasked050.ok
Resumable: a genome with both <gid>.assign.tsv and the sentinel is skipped; a genome that has
assign.tsv but no sentinel (i.e. produced before the second pass existed) gets the second pass
only, without re-running prodigal or the main scan.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_paths as P  # noqa: E402

SENTINEL = "unmasked050.ok"


def build_unmasked_db() -> Path:
    """Concatenate the group_4915 / group_5101 model.hmm files into one scan file.
    hmmsearch reads a plain multi-model HMM file directly; only hmmscan needs hmmpress."""
    out = P.CDS_DB_UNMASKED
    out.parent.mkdir(parents=True, exist_ok=True)
    srcs = [P.cds_profile(f) / "model.hmm" for f in P.UNMASKED_AT_050]
    for s in srcs:
        if not s.exists():
            raise SystemExit(f"missing profile HMM: {s}")
    if out.exists() and out.stat().st_mtime >= max(s.stat().st_mtime for s in srcs):
        return out
    with open(out, "wb") as fh:
        for s in srcs:
            fh.write(s.read_bytes())
    return out


def _best_from_domtbl(path: Path, best: dict) -> None:
    """domtblout: c0 target(protein) c3 query(family HMM) c5 qlen c7 full score c15 hmm_from c16 hmm_to."""
    for ln in path.read_text().splitlines():
        if ln.startswith("#") or not ln.split():
            continue
        c = ln.split()
        prot, fam = c[0], c[3]
        try:
            score = float(c[7]); hlen = int(c[5]); hf = int(c[15]); ht = int(c[16])
        except (ValueError, IndexError):
            continue
        cur = best.get(fam)
        if cur is None or score > cur[0]:
            best[fam] = (score, prot, hf, ht, hlen)


def _write_assign(assign: Path, best: dict) -> None:
    with open(assign, "w") as fh:
        fh.write("family\tprotein_id\thmm_cov\tscore\n")
        for fam, (score, prot, hf, ht, hlen) in best.items():
            cov = (ht - hf + 1) / hlen if hlen else 0
            fh.write(f"{fam}\t{prot}\t{cov:.3f}\t{score:.1f}\n")


def _read_assign(assign: Path) -> dict:
    """Re-read an existing assign.tsv as family -> (protein_id, hmm_cov, score) with the numeric
    fields kept as the original formatted strings, so _merge_rows() rewrites untouched families
    byte-for-byte."""
    rows = {}
    for i, ln in enumerate(assign.read_text().splitlines()):
        if i == 0 or not ln.strip():
            continue
        c = ln.split("\t")
        if len(c) >= 4:
            rows[c[0]] = (c[1], c[2], c[3])
    return rows


def _merge_rows(assign: Path, extra_best: dict) -> None:
    """Append/overwrite the extra families into an existing assign.tsv, preserving all other rows
    byte-for-byte."""
    rows = _read_assign(assign)
    for fam, (score, prot, hf, ht, hlen) in extra_best.items():
        cov = (ht - hf + 1) / hlen if hlen else 0
        rows[fam] = (prot, f"{cov:.3f}", f"{score:.1f}")
    with open(assign, "w") as fh:
        fh.write("family\tprotein_id\thmm_cov\tscore\n")
        for fam, (prot, cov, score) in rows.items():
            fh.write(f"{fam}\t{prot}\t{cov}\t{score}\n")


def _hmmsearch(db: str, faa: Path, dom: Path, env: dict) -> int:
    r = subprocess.run([P.HMMSEARCH, "--domtblout", str(dom), "--noali",
                        "-E", P.EVALUE, "--cpu", "2", db, str(faa)],
                       capture_output=True, text=True, env=env)
    return r.returncode


def predict_one(job):
    gid, path = job
    t0 = time.perf_counter()
    out = P.PREDICT / gid
    assign = out / f"{gid}.assign.tsv"
    sentinel = out / f"{gid}.{SENTINEL}"
    faa, ffn = out / f"{gid}.faa", out / f"{gid}.ffn"
    env = P.hmm_env()

    if assign.exists() and sentinel.exists():
        return (gid, "cached", 0.0, sum(1 for _ in open(assign)) - 1)

    if assign.exists() and faa.exists():
        # produced before the second pass existed: run only the two-model scan and merge
        dom2 = out / "dom_unmasked.txt"
        if _hmmsearch(str(P.CDS_DB_UNMASKED), faa, dom2, env) != 0:
            return (gid, "hmmsearch_unmasked_fail", time.perf_counter() - t0, 0)
        extra = {}
        _best_from_domtbl(dom2, extra)
        dom2.unlink(missing_ok=True)
        _merge_rows(assign, extra)
        sentinel.write_text("ok\n")
        return (gid, "topup", time.perf_counter() - t0, sum(1 for _ in open(assign)) - 1)

    out.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([P.PRODIGAL, "-i", path, "-a", str(faa), "-d", str(ffn),
                        "-p", "single", "-q", "-o", os.devnull, "-f", "gff"],
                       capture_output=True, text=True)
    if r.returncode != 0 or not faa.exists():
        return (gid, "prodigal_fail", time.perf_counter() - t0, 0)

    best = {}
    dom = out / "dom.txt"
    if _hmmsearch(str(P.CDS_DB), faa, dom, env) != 0:
        return (gid, "hmmsearch_fail", time.perf_counter() - t0, 0)
    _best_from_domtbl(dom, best)
    dom.unlink(missing_ok=True)

    dom2 = out / "dom_unmasked.txt"
    if _hmmsearch(str(P.CDS_DB_UNMASKED), faa, dom2, env) != 0:
        return (gid, "hmmsearch_unmasked_fail", time.perf_counter() - t0, 0)
    _best_from_domtbl(dom2, best)
    dom2.unlink(missing_ok=True)

    _write_assign(assign, best)
    sentinel.write_text("ok\n")
    return (gid, "ok", time.perf_counter() - t0, len(best))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(P.INPUT))
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    jobs = []
    with open(args.input) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi, pi = hdr.index("genome_id"), hdr.index("path")
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if len(c) > max(gi, pi):
                jobs.append((c[gi], c[pi]))
    if args.limit:
        jobs = jobs[:args.limit]

    P.PREDICT.mkdir(parents=True, exist_ok=True)
    build_unmasked_db()
    print(f"predict+assign over {len(jobs)} genomes ({args.workers} workers)", flush=True)
    t0 = time.perf_counter()
    ok = cached = topup = fail = done = 0
    hits, secs = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for gid, status, el, n in ex.map(predict_one, jobs):
            done += 1
            ok += status == "ok"; cached += status == "cached"; topup += status == "topup"
            if status in ("ok", "topup"):
                hits.append(n); secs.append(el)
            if status not in ("ok", "cached", "topup"):
                fail += 1
                sys.stderr.write(f"  {gid}: {status}\n")
            if done % 100 == 0:
                print(f"  {done}/{len(jobs)} (ok={ok} cached={cached} topup={topup} fail={fail})",
                      flush=True)
    wall = time.perf_counter() - t0
    print(f"done: {ok} predicted, {topup} topped up, {cached} cached, {fail} failed -> {P.PREDICT}")
    if hits:
        h = sorted(hits)
        n = len(h)
        print(f"HMM hits/genome: min={h[0]} p25={h[n//4]} median={h[n//2]} p75={h[3*n//4]} max={h[-1]}")
        s = sorted(secs)
        print(f"per-genome seconds (worker-side): median={s[len(s)//2]:.2f} max={s[-1]:.2f}")
    print(f"wall clock {wall:.1f}s for {len(jobs)} genomes = {wall/max(len(jobs),1):.3f} s/genome")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
