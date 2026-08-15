"""auric-gwas — command-line entry point.

Subcommands:
  scan      run the lineage-blind + lineage-aware association on a phenotype table (FUNCTIONAL)
  annotate  join variant_ids to the AURIC annotation layer (FUNCTIONAL)
  genotype  project user assemblies onto the frozen coordinate system (STUB -> cohort-genotyping fork)
  qc        filter user assemblies on the AURIC QC thresholds (STUB -> chop 03_qc.py)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import paths

paths.set_thread_caps(4)


def _cmd_scan(a: argparse.Namespace) -> int:
    from . import scan
    covar_cols = [c for c in (a.covar_cols.split(",") if a.covar_cols else []) if c]
    pheno_spec = dict(
        pheno_path=a.pheno_table, id_col=a.id_col, pheno_col=a.pheno_col,
        pheno_type=a.pheno_type, covar_cols=covar_cols,
        r_ge=a.r_ge, s_le=a.s_le, filter_col=a.filter_col, filter_val=a.filter_val,
    )
    out = scan.run(pheno_spec, source=a.source, matrix_dir=Path(a.matrix),
                   min_mac=a.min_mac, max_miss=a.max_miss, annotate_top=a.annotate_top)
    outdir = paths.assert_not_panel(a.out)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    res, diag = out["results"], out["diagnostics"]
    res.to_parquet(Path(outdir) / "association.parquet", index=False)
    res.head(a.annotate_top).to_csv(Path(outdir) / "association_top.tsv", sep="\t", index=False)
    (Path(outdir) / "diagnostics.json").write_text(json.dumps(diag, indent=2, default=float))
    print(json.dumps(diag, indent=2, default=float))
    print(f"\nwrote {outdir}/association.parquet ({len(res):,} variants), association_top.tsv, "
          f"diagnostics.json", file=sys.stderr)
    for w in diag["warnings"]:
        print(f"WARNING: {w}", file=sys.stderr)
    return 0


def _cmd_annotate(a: argparse.Namespace) -> int:
    from . import annotate
    ids = [l.strip() for l in Path(a.ids).read_text().splitlines() if l.strip()] if a.ids else a.variant_id
    df = annotate.annotate(ids)
    if a.out:
        df.to_csv(paths.assert_not_panel(a.out), sep="\t", index=False)
    else:
        print(df.to_string())
    return 0


def _cmd_genotype(a: argparse.Namespace) -> int:
    from . import genotype
    mdir = genotype.run_genotype(a.out, assemblies_dir=a.assemblies, manifest=a.manifest,
                                 workers=a.workers)
    print(f"cohort matrix: {mdir}\n"
          f"next: auric-gwas scan --source cohort --matrix {mdir} --pheno-table ...", file=sys.stderr)
    return 0


def _cmd_stub(name: str, pointer: str):
    def run(a):
        print(f"`{name}` is not yet implemented in the prototype.\n{pointer}", file=sys.stderr)
        return 2
    return run


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="auric-gwas", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="store_true", help="print version and the catalogue db_version")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("scan", help="lineage-blind + lineage-aware association")
    s.add_argument("--pheno-table", required=True, dest="pheno_table")
    s.add_argument("--id-col", required=True, dest="id_col", help="genome_id column")
    s.add_argument("--pheno-col", required=True, dest="pheno_col")
    s.add_argument("--pheno-type", choices=["binary", "continuous"], required=True, dest="pheno_type")
    s.add_argument("--r-ge", type=float, default=None, dest="r_ge", help="MIC >= this -> resistant (binary)")
    s.add_argument("--s-le", type=float, default=None, dest="s_le", help="MIC <= this -> susceptible (binary)")
    s.add_argument("--covar-cols", default=None, dest="covar_cols", help="comma-separated covariate columns")
    s.add_argument("--filter-col", default=None, dest="filter_col", help="subset the table, e.g. drug")
    s.add_argument("--filter-val", default=None, dest="filter_val")
    s.add_argument("--source", choices=["panel", "cohort"], default="panel")
    s.add_argument("--matrix", default=str(paths.MATRIX), help="genotype-matrix dir (frozen panel or a cohort build)")
    s.add_argument("--min-mac", type=int, default=5, dest="min_mac")
    s.add_argument("--max-miss", type=float, default=0.10, dest="max_miss")
    s.add_argument("--annotate-top", type=int, default=200, dest="annotate_top")
    s.add_argument("--out", required=True)
    s.set_defaults(func=_cmd_scan)

    an = sub.add_parser("annotate", help="annotate variant_ids against the frozen catalogue")
    g = an.add_mutually_exclusive_group(required=True)
    g.add_argument("--ids", help="file of variant_ids, one per line")
    g.add_argument("--variant-id", nargs="+", dest="variant_id")
    an.add_argument("--out", help="TSV output (default: stdout)")
    an.set_defaults(func=_cmd_annotate)

    gt = sub.add_parser("genotype", help="project assemblies onto the frozen coordinate system")
    src = gt.add_mutually_exclusive_group(required=True)
    src.add_argument("--assemblies", help="directory of assembly FASTAs (genome_id = file stem)")
    src.add_argument("--manifest", help="TSV with genome_id + path (or assembly_path) columns")
    gt.add_argument("--out", required=True, help="per-request output root (NOT under the frozen panel)")
    gt.add_argument("--workers", type=int, default=16)
    gt.set_defaults(func=_cmd_genotype)

    qc = sub.add_parser("qc", help="[stub] filter assemblies on the AURIC QC thresholds")
    qc.set_defaults(func=_cmd_stub(
        "qc",
        "Use scripts/chop_nicu_reanalysis/03_qc.py with config/config.yaml['qc'] thresholds "
        "(no dereplication for external cohorts)."))
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "version", False):
        from . import __version__
        print(f"auric-gwas {__version__}; catalogue db_version {paths.db_version()}")
        return 0
    if not getattr(args, "cmd", None):
        build_parser().print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
