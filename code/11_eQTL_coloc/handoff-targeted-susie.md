# Handoff: targeted `coloc.susie` analysis

## Objective

Design and run a maximally defensible targeted `coloc.susie` analysis on
`rsrv16`, beginning with the MAPK3 locus as a correctness and runtime pilot.
Run `susieR` directly as an independent implementation check. After the pilot,
prepare a concrete execution plan for the remaining targeted regions selected
from the approved final colocalization results.

This is exploratory work. Do not alter the manuscript tables or replace the
production `coloc.abf` analysis unless the investigators later make that
decision explicitly.

## Mandatory first step

Before executing the pilot, critically review this proposed workflow and
research the current official documentation. Correct any weak assumptions,
especially those involving LD, allele alignment, case-control summary data,
SuSiE prior settings, covariate adjustment, and summary-statistic/LD mismatch.
Explain any material correction before relying on it.

Primary documentation:

- [`coloc` SuSiE vignette](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html)
- [`coloc::runsusie()`](https://chr1swallace.github.io/coloc/reference/runsusie.html)
- [`coloc::coloc.susie()`](https://chr1swallace.github.io/coloc/reference/coloc.susie.html)
- [`susieR::susie_rss()`](https://stephenslab.github.io/susieR/reference/susie_rss.html)
- [`susieR` fine-mapping with summary statistics](https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html)
- [International Genome Sample Resource data](https://www.internationalgenome.org/data/)
- [1000 Genomes 30x high-coverage collection](https://internationalgenome.org/data-portal/data-collections/1000genomes_30x/)

Prefer official documentation, primary papers, and primary data sources. Record
the URLs, access date, package versions, and conclusions that affect the run.

## Current verified state

- Repository root on the server:
  `/home/gpertea/work/R/spatialDLPFC_mdd_bpd`
- Connect from home with `ssh rsrv16`. `srv16` is the direct/on-site alias.
- Branch: `devel`
- Verified HEAD: `45b5cfb4949ec8f6fae0a7ed62dceeb3e8b3dbfb`
- The server has no tracked Git modifications at handoff time.
- No SuSiE/coloc exploratory process is running.
- Server resources observed earlier: 32 CPU cores and approximately 503 GiB
  RAM. Neither `coloc.susie` nor `susieR` has a supported GPU backend; confirm
  this against current documentation and plan CPU parallelism across targets.

An isolated exploratory R library already exists at:

```text
processed-data/11_eQTL_coloc/seurat/coloc/susie_exploratory/Rlib
```

It contains official GitHub development versions installed without changing
the production R library:

- `coloc` 6.0.1, commit
  `50fe5291fea7f8ab49823bd86747385d6e56870f`
- `susieR` 0.16.6, commit
  `ef213feed2cb82419677661a8c986e1504df2c73`

The normal server library containing packages such as `pgenlibr` is:

```text
/scratch/gpertea/R/x86_64-pc-linux-gnu-library/4.5
```

Use both without upgrading the production library, for example:

```bash
export R_LIBS_USER="processed-data/11_eQTL_coloc/seurat/coloc/susie_exploratory/Rlib:/scratch/gpertea/R/x86_64-pc-linux-gnu-library/4.5"
```

Reconfirm the loaded package paths, versions, and Git commits in the run
metadata.

## Absolute GWAS data rule

Do **not** use the full-23andMe BD or MDD data. Those data were run only on a
separate branch and are not approved for this analysis.

Do not reuse the existing untracked `coloc/BD` or `coloc/MDD` caches or their
run metadata. Some old server-side caches refer to full-23andMe inputs. Their
presence does not make them valid.

Resolve the approved GWAS BCFs from the current `devel` version of
`code/11_eQTL_coloc/utils.R`. At handoff time the approved files are:

```text
BD/bip2024_eur_no23andMe.hg38.bcf
MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf
SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf
```

Verify the resolved absolute paths and checksums before analysis. Create fresh,
isolated targeted GWAS slices from those BCFs. Never write them into or infer
them from the old production cache directories.

## Inputs and target definition

Use the current approved final workbook only to define targets:

```text
processed-data/11_eQTL_coloc/seurat/final/coloc_results.xlsx
```

The `coloc_pass` sheet currently contains 100 targeted
disorder/context/gene comparisons. It includes five MAPK3 rows:

- BD: Astro, Inhb, L2.3, and L5
- MDD: Inhb

SCZ-MAPK3 is scientifically important as a distinct-signal/H3 comparison even
though it did not pass the final strong-H4 table. Include SCZ-MAPK3 in the
pilot. A practical runtime pilot may run the entire MAPK3 locus for all eight
all-donor eQTL contexts and all three disorders, while clearly separating the
five approved `coloc_pass` targets from diagnostic comparisons.

For the later all-target plan:

1. Start from the 100 rows in `coloc_pass`.
2. Remove MAPK3 comparisons already completed in the pilot.
3. Collapse overlapping gene windows into reusable genomic regions for input
   extraction and LD calculation, without collapsing distinct
   disorder/context/gene results.
4. Add H3/distinct-signal follow-ups only as a separately justified scientific
   target list. Do not silently expand to every tested gene.
5. Report counts of comparison rows, unique genes, unique regions, eQTL fits,
   GWAS fits, and reusable LD matrices before launching the full targeted run.

Do not launch the full all-target run until the MAPK3 validation and projected
resource plan have been reviewed with the user.

## Dense summary statistics

`coloc.susie` and `susie_rss` require dense regional summary statistics. Do
not use a GWAS-significance filter or the BH05-filtered eQTL table.

Use the full raw tensorQTL nominal parquet data on the server:

```text
processed-data/11_eQTL_coloc/seurat/tqtl_out/
```

Relevant files follow this pattern:

```text
<dataset>.gene.cis_qtl_pairs.chr<chromosome>.parquet
```

The all-donor contexts are `astro`, `inhb`, `l2-3`, `l4`, `l5`, `l6`,
`uvasc`, and `oligo`. Sample counts and covariate files are in:

```text
processed-data/11_eQTL_coloc/seurat/tqtl_in/
```

The established cis window is 1 Mb from the gene start/TSS. Recheck this in
the current code and preserve the exact tensorQTL variant universe. Read only
the required chromosome/gene rows; do not aggregate or copy the approximately
26 GiB parquet collection.

Create fresh dense GWAS slices at the exact candidate positions from the
approved no-23andMe/current SCZ BCFs. Preserve the current approved imputation
quality rule from the code, currently SI >= 0.8, after verifying field meaning
and availability for each disorder. Do not apply a GWAS p-value threshold.

## LD requirements

SuSiE requires a signed allele-aligned correlation matrix, `r`, not `r^2`.
Variant names and ordering must agree exactly between summary statistics and
the corresponding LD matrix.

### eQTL LD

The preferred source is the exact donor genotype data used by tensorQTL:

```text
processed-data/00_genotypes/plink2/merged_maf05.pgen
processed-data/00_genotypes/plink2/merged_maf05.pvar
processed-data/00_genotypes/plink2/merged_maf05.psam
```

Use the exact samples for each eQTL context. Most all-donor contexts have 119
samples; Oligo has fewer. Match sample IDs explicitly.

Critically determine from SuSiE theory and the tensorQTL model whether the LD
matrix should be computed from centered/scaled raw dosages or genotypes
residualized on the same context-specific covariate design used for the
summary-statistic regression. Do not assume either choice without documenting
the derivation. If uncertainty remains, run a declared sensitivity comparison.

### GWAS LD

Do not use the 119 brain donors as the GWAS LD reference. Prefer
study-matched LD if it can legitimately be obtained. Otherwise justify an
ancestry-matched external EUR reference for these EUR GWAS releases.

A targeted 1000 Genomes high-coverage GRCh38 EUR-unrelated MAPK3 extraction is
already present under:

```text
processed-data/11_eQTL_coloc/seurat/coloc/susie_exploratory/reference/
```

It contains 503 unrelated EUR samples and 1,941 biallelic SNVs overlapping the
MAPK3 tensorQTL positions. Treat it as a candidate input, not as automatically
validated truth. Recheck:

- sample selection and unrelated EUR status;
- intersection with the 3,202-sample high-coverage panel;
- GRCh38 coordinates;
- REF/ALT identity and ALT-dosage orientation;
- missingness, MAF, monomorphic variants, multiallelic handling, and checksums;
- whether 503 reference samples are adequate for each locus and allele-frequency
  range;
- summary-statistic/LD compatibility for each GWAS.

The high-coverage phased source file is:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

Use dosage correlations for SuSiE; phasing is not itself required for the
genotype-dosage correlation matrix. Do not conflate phased haplotype LD with
the signed dosage correlation required by `susie_rss`.

### LD and numerical diagnostics

At minimum, report:

- eigenvalue/positive-semidefinite checks and matrix symmetry;
- variant and sample counts before and after QC;
- allele and order concordance;
- `susieR::estimate_s_rss()` and appropriate `kriging_rss()` diagnostics;
- convergence, iteration count, credible-set purity, PIPs, and credible-set
  membership;
- missing or monomorphic variants and every exclusion reason.

Do not use `nearPD`, LD shrinkage, arbitrary pruning, or reduced credible-set
coverage merely to force a result. If exact or near-exact dosage ties cause a
numerical problem, preserve a reversible mapping from any computationally
collapsed representation back to **all** statistically indistinguishable
variant IDs. Never reintroduce the prior single-variant reporting loss.

## Required MAPK3 pilot

Use a fresh output directory, separate from both production outputs and the
aborted prototype, for example:

```text
processed-data/11_eQTL_coloc/seurat/coloc/susie_targeted_no23andMe/
```

The existing script
`susie_exploratory/04a_exploratory_MAPK3_coloc_susie.R` is an aborted prototype
and must not be executed as-is. It incorrectly selected GWAS caches from old
run metadata, thereby entering the full-23andMe path, and it failed during the
first ABF result extraction. Use it only as a warning or source of small,
independently reviewed utility fragments. This task is `coloc.susie` only; do
not rerun `coloc.abf`.

For every MAPK3 pilot comparison:

1. Build dense harmonized eQTL and approved GWAS datasets.
2. Build independently justified eQTL and GWAS signed LD matrices.
3. Run `susieR::susie_rss()` directly for each trait.
4. Run `coloc::runsusie()` on the exact same `z`, `R`, `N`, variant order, and
   SuSiE arguments.
5. Compare the direct `susieR` fit with the fine-mapping fit produced/annotated
   by `runsusie`: convergence, iterations, PIPs, credible sets, log Bayes
   factors, and variant ordering.
6. Pass the validated fits to `coloc::coloc.susie()` with the established
   cross-trait priors `p1 = 1e-4`, `p2 = 1e-4`, and `p12 = 1e-5`, unless a
   clearly labeled sensitivity analysis justifies alternatives.
7. Also verify that supplying prefit SuSiE objects versus a direct
   `coloc.susie()` call gives compatible results when the effective settings
   are identical.

`susieR` fine-mapping outputs and `coloc.susie` colocalization posteriors are
different estimands and should not be numerically identical. “Divergence” here
means the underlying SuSiE fits differ despite identical inputs/settings, or
the two supported `coloc.susie` routes disagree. If that occurs, stop scaling
and investigate argument translation, stale documentation, defaults, priors,
variant names/order, LD matrices, sample size, and package versions.

There is a known version-specific issue to resolve explicitly: the installed
`coloc::runsusie()` help mentions a `prior_variance` argument, while installed
`susieR` 0.16.6 exposes `scaled_prior_variance`. Do not pass an obsolete or
partially matched argument. Inspect current function signatures/source and
official documentation, establish the correct mapping, and record every
effective SuSiE parameter. Also determine rather than assume the correct
quantitative-trait scaling and case-control sample-size treatment.

Absence of a 95% credible set is a valid result. Diagnose it; do not lower
coverage simply to manufacture a colocalization result.

## Timing and scaling report

Separate and record time for:

- targeted parquet/BCF extraction;
- genotype loading and LD construction;
- each direct `susie_rss` fit;
- each `runsusie` fit;
- each `coloc.susie` comparison;
- serialization and diagnostics.

Also record peak memory, number of variants, number of samples, number of
credible sets, convergence status, and failures. Use the MAPK3 measurements to
project total CPU time, wall time, storage, and safe parallelism for the
remaining target manifest. Cache and reuse only inputs/fits whose defining
disorder, context, region, variants, LD source, and settings are truly
identical.

## Reproducibility and stopping conditions

- Keep all exploratory outputs isolated from production coloc directories.
- Do not modify the final XLSX workbooks.
- Do not run tensorQTL, regenerate eQTL tables, rerun `coloc.abf`, or alter the
  newly updated DEG annotations.
- Write a target manifest, input provenance table, variant-exclusion audit,
  parameter table, timing table, diagnostics table, session information, and
  checksums.
- Checkpoint per comparison so a partial failure is recoverable.
- If LD provenance or allele harmonization cannot be made defensible, stop and
  report the specific blocker rather than producing a posterior from uncertain
  inputs.
- After the MAPK3 pilot, report findings and the all-target execution plan to
  the user. Wait for approval before launching the remaining regions.
