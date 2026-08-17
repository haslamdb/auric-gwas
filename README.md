# auric-gwas

This is an end to end workflow that projects *Staphylococcus aureus* genomes onto the frozen **AURIC** coordinate system and runs
**lineage-aware genome-wide association**. AURIC assigns every core single-nucleotide variant a stable
identifier defined by gene family, profile-HMM match state and codon position, so genomes added years
apart — or from different cohorts — share one coordinate system. `auric-gwas` genotypes your assemblies
onto that system (read-only against the frozen reference) and runs the two-pass association the AURIC
paper recommends: a lineage-unadjusted scan and a lineage-aware mixed model (EMMAX-style, with an
AURIC-derived kinship matrix), returning every variant ranked and annotated.

## Why

Bacterial GWAS is dominated by clonal population structure: resistance and thousands of unrelated
alleles are inherited together within lineages, so an unadjusted scan reports large numbers of variants
that merely mark the resistant clones. The AURIC-derived kinship is computed once from the frozen
catalogue and released with it, so a variant's lineage-adjusted significance means the same thing in
every study that uses it — a standardized correction for clonal confounding.

## Install

```bash
conda env create -f environment.yml     # analysis deps + the genotyping tool chain on PATH
conda activate auric-gwas
pip install -e .
```

`auric-gwas` needs the **AURIC reference bundle** (the frozen profiles, coordinate table and catalogue;
distributed separately — a Zenodo deposit aligned with the manuscript). Unpack it and point `AURIC_HOME`
at it (default `/fastpool/sausnp` on the development machine):

```bash
export AURIC_HOME=/path/to/auric_reference
auric-gwas --version        # prints the tool version and the catalogue db_version
```

## Pipeline

```bash
# 1. Project assemblies onto AURIC's 1,231,021-column coordinate system (read-only; the frozen
#    reference is never modified). Produces a cohort genotype matrix in the frozen column order.
auric-gwas genotype --assemblies my_isolates/ --out run/

# 2. Lineage-unadjusted + lineage-aware association against a phenotype, with a covariate.
auric-gwas scan \
    --source cohort --matrix run/genotype_matrix \
    --pheno-table phenotypes.tsv --id-col genome_id \
    --pheno-col mic --pheno-type binary --r-ge 4 --s-le 1 \
    --covar-cols cohort --out results/
```

`scan` writes `association.parquet` (every variant, lineage-blind and lineage-aware ranks, joined to
AURIC annotation), `association_top.tsv`, and `diagnostics.json` (λ_GC for both scans, h², the
significant-count collapse, and over-correction warnings). A continuous phenotype (e.g. log₂ MIC) uses
`--pheno-type continuous`; covariates are `--covar-cols a,b`.

## Worked example (reproduces the paper's positive control)

Levofloxacin, 266 clinical-MIC isolates, cohort as a fixed covariate. The scan is **blind to functional
annotation** — no locus is named in advance — yet recovers the quinolone determinant *gyrA* Ser84Leu as
the top association:

```
samples: 266  (118 R / 148 S)
lineage-blind significant: 16,205   ->   lineage-aware significant: 10
genomic-control lambda: blind 22.1  aware 0.010   (h2_kinship 0.998)

top lineage-aware associations:
 rank_aware gene_symbol aa_change     effect   p_aware    p_blind    n_origins
          1        gyrA      S84L   missense  5.7e-32    3.6e-129        199
```

`grlA`/`parC` Ser80Tyr follows at rank 7. At this sample size the phenotype is nearly perfectly
lineage-confounded (h²=0.998), so the mixed model over-corrects (λ_aware=0.010) and the tool says so —
the aware ranking is conservative, and coding consequence + genome context (both in the annotated
output) separate a determinant from a lineage marker.

Reproduce it (needs `AURIC_HOME` and the `mic_phenotypes.tsv` example phenotype):

```bash
AURIC_HOME=/path/to/auric_reference python auric_gwas/examples/levofloxacin.py
```

## Subcommands

| command | status | does |
|---|---|---|
| `genotype` | ✓ | project assemblies onto the frozen coordinate system → cohort matrix |
| `scan` | ✓ | lineage-blind + lineage-aware association on a phenotype table |
| `annotate` | ✓ | join variant_ids to the AURIC annotation layer |
| `qc` | stub | filter assemblies on the AURIC QC thresholds (length, contigs, N50, CheckM2, ANI) |

## Read-only guarantee

`auric-gwas` never writes under the AURIC reference. The genotyping stages run through a guard that
asserts every write target is under your `--out` root and refuses to run if any resolves under the
reference (`$AURIC_HOME/{db,profiles,profiles_masked}`); the catalogue is opened `read_only=True`. The
projection is validated bit-identical (concordance 1.000000) against reference genomes.

## Method

The lineage-aware scan is the EMMAX approach (Kang et al., *Nat Genet* 2010): one variance-component
ratio is estimated from the kinship-only null by REML, then held fixed while every variant is tested by
GLS. The genotyping projects each assembly onto the frozen profile HMMs (prodigal → hmmsearch →
hmmalign → match-state→codon walk), so a genome subset reproduces exactly the alleles the full reference
run produced.

## License

MIT © 2026 David B. Haslam
