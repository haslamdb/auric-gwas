# auric-gwas

This is an end-to-end workflow that finds the genetic changes associated with a trait, including antibiotic
resistance and virulence across a collection of *Staphylococcus aureus* genomes, while
correcting for the fact that bacteria inherit most of their DNA in a block.

`auric-gwas` places your assembled genomes onto **AURIC**, a frozen, species-wide coordinate system
for *S. aureus*, and runs a **genome-wide association study (GWAS)**: it tests every catalogued
single-nucleotide variant for correlation with your phenotype, then re-tests it with a correction for
shared ancestry. The output is every variant, ranked twice and joined to AURIC's annotation, so you
can tell a genuine determinant from a passenger that merely marks a lineage.

---

## The problem, in plain terms

Two facts about bacteria make GWAS hard, and AURIC addresses both.

**1. Bacteria are clonal, so almost everything is correlated with almost everything.** *S. aureus*
reproduces by copying itself. A successful strain leaves thousands of near-identical descendants that
share the whole genome, not just the one mutation that matters. So a strain that is resistant to a
drug also carries thousands of *unrelated* mutations it happened to inherit from the same ancestor. A
naive association test flags all of them. In the worked example below, an uncorrected scan calls
**16,205** variants "significant" for fluoroquinolone resistance; almost all of them simply mark the
resistant lineages rather than cause resistance. This is **clonal (population-structure) confounding**,
and controlling it is the central task of any bacterial GWAS. The standard fix is a **kinship matrix**
(a table of how related every pair of strains is) fed to a **mixed model** that discounts similarity
two strains share only because they are cousins.

**2. There is no stable place to put a variant.** The usual way to name a variant is by its position
on one chosen **reference genome** — "position 1,473,201 of NCTC 8325." This breaks down for a
species like *S. aureus*, organized into ~55 divergent lineages (**clonal complexes**), for two
reasons. First, genes present in your strains but absent from the reference simply cannot be named at
all. Second, every time a study rebuilds its variant set from a fresh genome comparison, the positions
renumber — so variant "position 41,203" in one paper and the same biological change in another paper
have different names, and results cannot be pooled or compared across studies.

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
- **Honest positions only.** A position is scored only where **orthology** is guaranteed — the
  position is genuinely "the same position by descent" across genomes: present in nearly all strains,
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
| **Population-structure correction** | You provide it (a phylogeny or a distance/similarity matrix). | Computed for you as a genomic relatedness matrix from your cohort's genotypes *on the AURIC coordinates*, then applied EMMAX-style. |
| **Annotation** | Separate step you run afterward on your hits. | Top hits arrive pre-joined to AURIC's layer: coding consequence, amino-acid change, `n_origins`, ESM-2 constraint, AMR/virulence labels, recombination fraction. |
| **Determinant vs. lineage marker** | You infer it yourself. | The homoplasy count (`n_origins`) and coding consequence are in the output to make that call. |
| **Scope** | Any bacterium, any variant representation. | *S. aureus* only; core single-nucleotide variants only. |

The short version: use `pyseer` when you need k-mer/accessory-genome association, a non-model
organism, or full control over the variant representation. Use `auric-gwas` when your organism is
*S. aureus*, you want core-SNP association on a stable coordinate system, and you want results that
pool cleanly with other AURIC-based studies and come annotated for interpretation.

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
   matrix (GRM) is built from your cohort's genotypes, one variance-component ratio is estimated once
   from the kinship-only null by REML, then held fixed while every variant is tested by generalized
   least squares. This discounts the similarity strains share by descent.

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
  genomic-control **λ_GC** for both scans (how inflated each is), the kinship heritability **h²**, the
  significant-count collapse from blind to aware, and explicit **over-correction warnings**.

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
