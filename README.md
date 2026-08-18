# auric-gwas

An end-to-end workflow for associating genetic variants with a trait — antibiotic resistance,
virulence, host — across a collection of *Staphylococcus aureus* genomes, with the clonal population
structure that dominates bacterial genomes controlled by a lineage-aware mixed model.

`auric-gwas` places your assembled genomes onto **AURIC**, a frozen, species-wide coordinate system
for *S. aureus*, and runs a **genome-wide association study (GWAS)**: it tests every catalogued
single-nucleotide variant for correlation with your phenotype, then re-tests it with a correction for
shared ancestry. The output is every variant, ranked twice and joined to AURIC's annotation, so you
can tell a genuine determinant from a passenger that merely marks a lineage.

---

## Two problems specific to bacterial GWAS

**1. Clonal population structure inflates the association.** *S. aureus* is a largely clonal,
asexually reproducing species: it recombines rarely relative to its mutation rate, so linkage
disequilibrium extends across the entire chromosome rather than decaying over kilobases as it does in
sexually recombining organisms. A resistant strain and its recent descendants therefore share not
only the causal allele but the tens of thousands of neutral, co-inherited variants that mark their
lineage — the genomic background and the phenotype are confounded by shared ancestry. An unadjusted
scan cannot separate the two: in the worked example below it returns **16,205** genome-wide-significant
variants for fluoroquinolone resistance, nearly all of them lineage markers rather than determinants.
The correction is a **mixed model** with a **kinship (genomic relatedness) matrix** as a random
effect, which absorbs the covariance two isolates share by descent and tests each variant against that
background.

**2. Variant coordinates are not stable across studies.** The conventional identifier is a position on
one chosen **reference genome** — "position 1,473,201 of NCTC 8325." For a species partitioned into
~55 divergent lineages (**clonal complexes**) with a large accessory genome, that fails twice. Core
genes absent from the reference cannot be addressed at all, and every study that recalls variants
against a different reference — or rebuilds its coordinate set from a fresh pangenome — renumbers
everything, so the same biological substitution carries a different identifier in each paper and
results cannot be pooled or meta-analyzed.

---

## What AURIC is, and why it helps

**AURIC** (**AU**reus **R**eference-**I**nvariant **C**oordinates) is a frozen reference catalogue of
single-nucleotide variation in *S. aureus*, built in the companion project
[`staph_snp_db`](../staph_snp_db) from **20,685 quality-controlled genomes** (4,231 distinct
lineages). It solves the naming problem by numbering positions on **statistical templates of each
gene** instead of on one chromosome:

- Each core gene is turned into a **profile hidden Markov model (profile HMM)** — a consensus template
  whose numbered columns are called **match states**.
- A variant's address is therefore `(gene family, match state, codon position)` — for example "*gyrA*,
  match state 84, first codon base" — not a coordinate on a single genome.
- Because a new genome is *compared* to the templates rather than used to rebuild them, adding genomes
  later never moves a column. The templates were **frozen** once built; every variant identifier
  (`SAU_…`) is permanent and is never renumbered.

That frozen coordinate system is what makes AURIC useful beyond one study:

- **Portability.** A variant means the same thing in your cohort, in a collaborator's cohort, and in
  the catalogue. Genomes sequenced years apart, or by different labs, share one coordinate system, so
  results can be pooled and compared directly.
- **Completeness.** The catalogue spans the whole species' core genome (1,231,021 high-confidence
  variant positions), including genes absent from any single reference strain.
- **Orthology-guaranteed positions only.** A position is scored only where **orthology** holds — where
  it is genuinely "the same position by descent" across genomes: present in nearly all strains,
  single-copy, and not inside a repeat. Everything else is **masked** and reported, not silently
  dropped.
- **A homoplasy statistic.** For every variant AURIC records **`n_origins`** — how many times the
  change arose *independently* across the species tree. A mutation that arose once and rode a lineage's
  expansion (`n_origins = 1`) is uninformative; one that arose independently 199 times (like *gyrA*
  Ser84Leu, below) is a strong signature of selection. This is exactly the axis that separates a
  resistance determinant from a lineage marker, and `auric-gwas` puts it in the output.

For the full build story, written for non-specialists, see
[`staph_snp_db/docs/AURIC_background_and_build.md`](../staph_snp_db/docs/AURIC_background_and_build.md).

---

## How this differs from pyseer and other bacterial GWAS tools

`pyseer` (and `treeWAS`, `Scoary`, `bugwas`, `hogwash`) are general bacterial-GWAS engines: you bring
your own variants — k-mers, unitigs, or SNPs you called against a reference you chose — and your own
population-structure estimate. They are flexible and species-agnostic. `auric-gwas` is narrower and
species-specific, and trades that flexibility for standardization.

| | pyseer / general tools | auric-gwas |
|---|---|---|
| **Variant space** | Whatever you supply — k-mers, unitigs, or SNPs vs. your chosen reference. Identifiers are study-local. | The frozen AURIC catalogue. Identifiers (`SAU_…`) are permanent and shared across every study. |
| **Cross-study comparability** | Limited — a k-mer or a reference position means different things in different runs. | Built in — the coordinate system is the same catalogue for everyone. |
| **Population-structure correction** | You supply the similarity matrix (from a phylogeny or a core-genome distance) and choose a model — fixed-effects with MDS covariates, a FaST-LMM mixed model (`--lmm`), or elastic net. | Fixed choice: a genomic relatedness matrix built for you from your cohort's genotypes *on the AURIC coordinates*, applied EMMAX-style. |
| **Annotation** | Separate step you run afterward on your hits. | Top hits arrive pre-joined to AURIC's layer: coding consequence, amino-acid change, `n_origins`, ESM-2 constraint, AMR/virulence labels, recombination fraction. |
| **Determinant vs. lineage marker** | You infer it yourself. | The homoplasy count (`n_origins`) and coding consequence are in the output to make that call. |
| **Scope** | Any bacterium, any variant representation. | *S. aureus* only; core single-nucleotide variants only. |

**On the mixed model specifically.** When pyseer is run in its LMM mode (`--lmm`), the statistical
model is essentially the one `auric-gwas` uses. Both fit a single-kinship linear mixed model
(y = Xβ + u + ε, u ~ N(0, σ²_g·K)) and both estimate the variance component **once** under the null
— covariates and kinship, no tested variant — then hold it fixed while scanning the genome, rather
than re-fitting per variant. They differ in lineage and detail, not in kind: EMMAX (Kang et al.,
*Nat Genet* 2010) eigendecomposes K, estimates the variance ratio by REML, and applies a rotated GLS
Wald test per variant; pyseer's LMM is FaST-LMM (Lippert et al., *Nat Methods* 2011), which estimates
the component by ML and tests with a likelihood-ratio statistic, its low-rank factorization giving a
speed advantage when few features build K. On a full-rank GRM the two produce nearly identical
p-values. The consequential difference between the tools is therefore the **input** — pyseer LMM
typically over k-mers/unitigs with a phylogeny-derived similarity matrix, `auric-gwas` over AURIC
core-SNP genotypes with a GRM built from them — not the mixed-model math. pyseer's other modes
(fixed-effects with MDS covariates, elastic net) are genuinely different approaches; EMMAX
corresponds to `--lmm`, not to those.

The short version: use `pyseer` when you need k-mer/accessory-genome association, a non-model
organism, or full control over the variant representation and the association model. Use `auric-gwas`
when your organism is *S. aureus*, you want core-SNP association on a stable coordinate system, and
you want results that pool cleanly with other AURIC-based studies and come annotated for
interpretation.

---

## Install

```bash
conda env create -f environment.yml     # analysis deps + the genotyping tool chain on PATH
conda activate auric-gwas
pip install -e .
```

`auric-gwas` needs the **AURIC reference bundle** (the frozen profiles, coordinate table, and
catalogue; distributed separately as a Zenodo deposit aligned with the manuscript). Unpack it and
point `AURIC_HOME` at it (default `/fastpool/sausnp` on the development machine):

```bash
export AURIC_HOME=/path/to/auric_reference
auric-gwas --version        # prints the tool version and the catalogue db_version
```

---

## The workflow

`auric-gwas` runs in two commands: **genotype** your assemblies onto the coordinate system, then
**scan** the resulting matrix against a phenotype.

### Step 1 — `genotype`: place your genomes on the coordinate system

```bash
auric-gwas genotype --assemblies my_isolates/ --out run/
```

Each assembly is read gene by gene: open reading frames are predicted (**Prodigal**), matched to the
frozen profile HMMs (**hmmsearch**), aligned to their template (**hmmalign**), and each nucleotide is
walked from match state to codon position. The result is a **cohort genotype matrix** — one row per
genome, one column per AURIC variant, in the frozen column order — written under `run/`. Because the
templates are never modified, a genome scored in isolation reproduces exactly the allele it would get
in the full 20,685-genome reference run (validated bit-identical, concordance 1.000000).

### Step 2 — `scan`: association, run twice

```bash
auric-gwas scan \
    --source cohort --matrix run/genotype_matrix \
    --pheno-table phenotypes.tsv --id-col genome_id \
    --pheno-col mic --pheno-type binary --r-ge 4 --s-le 1 \
    --covar-cols cohort --out results/
```

`scan` tests every variant twice against your phenotype:

1. **Lineage-blind** — an ordinary linear model, no ancestry correction. This is the scan dominated by
   clonal confounding; it is reported so you can *see* the confounding, not because you should trust it.
2. **Lineage-aware** — an **EMMAX** mixed model (Kang et al., *Nat Genet* 2010). A genomic relatedness
   matrix (GRM) is built from your cohort's genotypes, and one variance-component ratio (the balance of
   lineage variance to residual variance) is estimated once from the kinship-only null by **REML**
   (restricted maximum likelihood — the standard estimator for variance components in a mixed model,
   which corrects the downward bias of ordinary maximum likelihood). That ratio is then held fixed
   while every variant is tested by **GLS** (generalized least squares — least squares that weights
   observations by the inverse of their covariance, so correlated isolates are down-weighted). This
   discounts the similarity strains share by descent.

Phenotypes can be **binary** (e.g. resistant vs. susceptible from an MIC breakpoint via `--r-ge` /
`--s-le`) or **continuous** (`--pheno-type continuous`, e.g. log₂ MIC). Covariates are
`--covar-cols a,b` (cohort, sequencing platform, collection year, …).

### Outputs

`scan` writes to `results/`:

- **`association.parquet`** — every tested variant, its blind and aware p-values and ranks, the
  aware effect size, joined to AURIC annotation (gene, amino-acid change, `n_origins`, consequence,
  AMR labels).
- **`association_top.tsv`** — the top hits, for quick reading.
- **`diagnostics.json`** — the health checks that tell you whether to trust the run: the
  **genomic-control λ_GC** for both scans, the kinship heritability **h²**, the significant-count
  collapse from blind to aware, and explicit **over-correction warnings**. λ_GC is the ratio of the
  observed median test statistic to the value expected under the null (Devlin & Roeder, 1999); λ_GC ≈ 1
  is calibrated, λ_GC > 1 indicates inflation from uncontrolled population structure, and λ_GC < 1
  indicates over-correction. h² is the proportion of phenotypic variance explained by the kinship
  matrix — the fraction attributable to lineage; h² near 1 means the phenotype is almost entirely
  predicted by ancestry, so the mixed model has little residual signal to work with and tends to
  over-correct.

Read `diagnostics.json` before the hit list. A λ_aware far below 1 means the mixed model
over-corrected (see the example), in which case the aware ranking is conservative and you lean on
`n_origins` and coding consequence to interpret.

---

## Practical example (reproduces the paper's positive control)

Levofloxacin resistance, 266 clinical-MIC isolates, cohort as a fixed covariate. The scan is **blind
to functional annotation** — no locus is named in advance — yet recovers the known fluoroquinolone
determinant *gyrA* Ser84Leu as the top association:

```
samples: 266  (118 R / 148 S)
lineage-blind significant: 16,205   ->   lineage-aware significant: 10
genomic-control lambda: blind 22.1  aware 0.010   (h2_kinship 0.998)

top lineage-aware associations:
 rank_aware gene_symbol aa_change     effect   p_aware    p_blind    n_origins
          1        gyrA      S84L   missense  5.7e-32    3.6e-129        199
```

The blind scan's λ_GC of **22.1** is the confounding made visible: it flags 16,205 variants. The aware
scan collapses that to **10** and puts the true determinant first — a variant that arose independently
**199 times** (`n_origins = 199`), the homoplasy signature of selection. `grlA`/`parC` Ser80Tyr, the
other classic quinolone determinant, follows at rank 7.

At this sample size the phenotype is nearly perfectly lineage-confounded (h² = 0.998), so the mixed
model **over-corrects** (λ_aware = 0.010) and the tool says so in `diagnostics.json`. When that
happens the aware ranking is conservative rather than final, and the coding consequence plus `n_origins`
(both in the annotated output) are what separate a determinant from a lineage marker.

Reproduce it (needs `AURIC_HOME` and the `mic_phenotypes.tsv` example phenotype):

```bash
AURIC_HOME=/path/to/auric_reference python auric_gwas/examples/levofloxacin.py
```

---

## Subcommands

| command | status | does |
|---|---|---|
| `genotype` | ✓ | project assemblies onto the frozen coordinate system → cohort matrix |
| `scan` | ✓ | lineage-blind + lineage-aware association on a phenotype table |
| `annotate` | ✓ | join variant_ids to the AURIC annotation layer |
| `qc` | stub | filter assemblies on the AURIC QC thresholds (length, contigs, N50, CheckM2, ANI) |

---

## Read-only guarantee

`auric-gwas` never writes under the AURIC reference. The genotyping stages run through a guard that
asserts every write target is under your `--out` root and refuses to run if any resolves under the
reference (`$AURIC_HOME/{db,profiles,profiles_masked}`); the catalogue is opened `read_only=True`. The
projection is validated bit-identical (concordance 1.000000) against reference genomes.

---

## Method

The lineage-aware scan is **EMMAX** (Kang et al., *Nat Genet* 2010): one variance-component ratio is
estimated from the kinship-only null by REML, then held fixed while every variant is tested by GLS —
the "eXpedited" shortcut that avoids re-fitting a full mixed model per variant. Significance uses a
pattern-collapsed Bonferroni threshold (distinct genotype patterns, matching pyseer's convention).
Genotyping projects each assembly onto the frozen profile HMMs (Prodigal → hmmsearch → hmmalign →
match-state→codon walk), so a genome subset reproduces exactly the alleles the full reference run
produced.

The reference catalogue itself was built and frozen in the companion project
[`staph_snp_db`](../staph_snp_db); `auric-gwas` is the downstream, read-only association tool.

## License

MIT © 2026 David B. Haslam
