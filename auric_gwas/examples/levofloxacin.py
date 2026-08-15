#!/usr/bin/env python3
"""Worked example: recover the fluoroquinolone determinant from an AURIC lineage-aware scan.

Reproduces the manuscript's positive control through the AURIC-associate API: a lineage-unadjusted
scan calls thousands of variants driven by clonal structure, while the lineage-aware mixed model
(AURIC-derived kinship + cohort covariate) collapses them to the QRDR determinants, gyrA Ser84Leu
first. This is the example figure/box for the paper — the tool finds the known site without being told
where it is.

Data: results/analysis/amr/mic_phenotypes.tsv (levofloxacin clinical MICs, provenance-checked).
Run:  AURIC_HOME=/fastpool/sausnp /home/david/miniforge3/bin/python auric_gwas/examples/levofloxacin.py
"""
import os
import sys

from auric_gwas import scan

# The levofloxacin MIC table ships with the AURIC reference bundle / manuscript supplement,
# not with the code. Pass its path as argv[1] or set AURIC_PHENO.
PHENO = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "AURIC_PHENO", os.path.join(os.environ.get("AURIC_HOME", "."), "examples", "mic_phenotypes.tsv"))

out = scan.run(
    dict(
        pheno_path=PHENO,
        id_col="genome_id", pheno_col="mic", pheno_type="binary",
        r_ge=4, s_le=1,                       # CLSI levofloxacin breakpoints
        covar_cols=["cohort"],
        filter_col="drug", filter_val="Levofloxacin",
    ),
    source="panel",                           # analyse catalogue genomes; use source="cohort" for your own
)
res, diag = out["results"], out["diagnostics"]

print(f"samples: {diag['n_samples']}  ({diag['n_case']} R / {diag['n_control']} S)")
print(f"lineage-blind significant: {diag['n_sig_blind']:,}   ->   lineage-aware significant: {diag['n_sig_aware']}")
print(f"genomic-control lambda: blind {diag['lambda_blind']:.1f}  aware {diag['lambda_aware']:.3f}"
      f"   (h2_kinship {diag['h2_kinship']:.3f})\n")

cols = ["rank_aware", "gene_symbol", "aa_change", "effect", "p_aware", "p_blind", "n_origins"]
print("top lineage-aware associations:")
print(res.sort_values("p_aware").head(5)[cols].to_string(index=False))

for w in diag["warnings"]:
    print(f"\nNOTE: {w}")
