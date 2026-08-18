# eQTL and Colocalization Workflow

This folder contains the MBv eQTL and colocalization workflow.

## Concurrent Work

Assume other users, agents, editors, and jobs may change files at any time.

Only touch files and metadata needed for the current request. Never clean up,
revert, normalize, or explain unrelated changes unless they block the task.

Default scope:

- input annotation: Seurat label-transfer clusters
- samples: all 119 donors together
- split group: `all` only
- tensorQTL series: `seurat`
- DEG support: `broad_interaction`

By default, eQTL and colocalization analyses produce tables only for the
all-donor `all` group using Seurat clusters. Additional models, custom clusters,
PRECAST details, sex-stratified eQTLs, and other DEG views are documented in
`README_extra_models.md`.

## Inputs and Outputs

Canonical inputs:

- `processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata`
- `processed-data/00_genotypes/plink2/merged_maf05.{pgen,psam,pvar}`
- `processed-data/00_genotypes/plink2/merged_maf05_pca.eigenvec`
- `processed-data/ref/granges.qs2`
- `processed-data/ref/GWAS/{BD,MDD,SCZD}`

Canonical Seurat eQTL paths:

- `processed-data/11_eQTL_coloc/seurat/tqtl_in`
- `processed-data/11_eQTL_coloc/seurat/tqtl_out`
- `processed-data/11_eQTL_coloc/seurat/tables`
- `processed-data/11_eQTL_coloc/seurat/coloc`
- `plots/11_eQTL_coloc/coloc`

Default Seurat clusters:

- `Astro`
- `Inhb`
- `L2.3`
- `L4`
- `L5`
- `L6`
- `Micro.Vasc` (`uvasc` in dataset IDs)
- `Oligo`

## Workflow Scripts

Major workflow files:

- `stage_required_data.sh`: stage required local inputs.
- `00_get_SNP_PCs.sh`: compute genotype PCs.
- `01_prep_inputs.R`: prepare all-donor Seurat tensorQTL inputs.
- `02_run_tensorQTL.sh`: run tensorQTL for prepared datasets.
- `02a_tensorQTL_cis.py`: tensorQTL cis wrapper called by the run script.
- `03_eqtl_explore.Rmd`: build primary cis, independent, and significant-pair eQTL summary tables.
- `03ab_recover_indepeQTLs.R`: recover every exact dosage-equivalent variant tied with a reported cis or conditionally independent eQTL signal and write the expanded signal-member tables used by downstream exact-variant analyses.
- `README_eQTL_exact_tie_recovery_impact.md`: document the recovery scale and its eGene, GWAS, trifecta, and strong-coloc effects by disorder and domain-CT.
- `03a_nominal_eQTLs.Rmd`: gather nominal all-donor Seurat eQTLs.
- `03b_eQTL_boxplots.Rmd`: select example eQTL pairs and render genotype boxplots.
- `04_run_coloc.R`: run coloc ABF and sensitivity checks from full tensorQTL nominal parquet.
- `05_coloc_explore.Rmd`: flatten coloc outputs, apply sensitivity gates, write tables, and plot strong coloc counts.
- `06_eqtl_coloc_final_tables.Rmd`: rebuild the committed manuscript-ready eQTL and coloc workbooks from the recovered exact-tie table and final coloc outputs.
- `check_datatable_scoping.R`: heuristic scan for risky bare-symbol data.table joins.
- `utils.R`: shared DEG, GWAS, tensorQTL, summary, and plotting helpers.

## R data.table Safety

Bare names inside `DT[...]` should mean data.table columns. Do not create
caller-scope variables whose names match any column in a data.table used in
that same scope. This is a general rule, not limited to eQTL key columns.

This failure mode is silent: `DT[col == col]`, `DT[col %in% col]`, and
`DT[J(col)]` can use the column twice when a caller scalar named `col` was
intended. Mixed-symbol filters are just as dangerous: `DT[col == var_id]` can
silently use column `var_id` when the coding user or agent intended a caller
variable named `var_id`. Numeric results may look plausible while labels,
filters, or joins are wrong.

Static scans cannot prove safety when table columns come from loaded files.
Runtime object introspection is required because only the live data.table knows
whether a bare RHS/filter symbol is also a column.


## Prepare Inputs

Stage required local inputs as needed if working locally instead of on the cluster:

```bash
./stage_required_data.sh
```

Compute genotype PCs:

```bash
./00_get_SNP_PCs.sh
```

Prepare Seurat all-donor tensorQTL inputs:

```bash
Rscript ./01_prep_inputs.R
```

Validation without writing outputs:

```bash
Rscript ./01_prep_inputs.R --check-only
```

Default covariates:

```text
DX + sex + age + PC3 + nspots + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + expression PCs
```

Prepared all-donor dataset IDs have no split suffix, for example `astro`, `inhb`, `l2-3`, and `uvasc`.

`01_prep_inputs.R` writes `processed-data/11_eQTL_coloc/seurat/tqtl_in/prep_manifest.csv`.
This manifest records dataset IDs, split labels, context labels, and `n_samples`.
Downstream coloc uses the manifest for eQTL sample size `N`; the coloc SNP-level
inputs still come from full nominal parquet files in `seurat/tqtl_out`.

## Run tensorQTL

Run all prepared Seurat all-donor datasets:

```bash
./02_run_tensorQTL.sh
```

List commands without running:

```bash
./02_run_tensorQTL.sh --dry-run
```

Run one dataset:

```bash
./02_run_tensorQTL.sh astro
```

Resume from a dataset:

```bash
./02_run_tensorQTL.sh --start-from l5
```

The wrapper reads `prep_manifest.csv` when present and defaults to `--splits all`.

tensorQTL nominal outputs for coloc are the per-chromosome parquet files:

```text
processed-data/11_eQTL_coloc/seurat/tqtl_out/<dataset_id>.gene.cis_qtl_pairs.chr*.parquet
```

They are generated by tensorQTL nominal mapping with statistics written for all
cis gene-SNP pairs, not only eGenes or significant pairs.

## Run Coloc ABF

Run MDD and BD coloc ABF for all Seurat all-donor domainCTs:

```bash
Rscript ./04_run_coloc.R
```

Run one disorder only, useful for separate machines:

```bash
Rscript ./04_run_coloc.R --disorder MDD
Rscript ./04_run_coloc.R --disorder BD
```

Optional filters:

```bash
Rscript ./04_run_coloc.R --datasets astro,l2-3 --disorder MDD
Rscript ./04_run_coloc.R --disorder BD --chromosomes chr22
Rscript ./04_run_coloc.R --dry-run --disorder BD
Rscript ./04_run_coloc.R --disorder MDD --n-cores 12
```

Coloc defaults:

- disorders: `MDD`, `BD`
- datasets: Seurat all-donor `split == "all"` domainCTs
- loci: all genes present in tensorQTL nominal parquet after GWAS/eQTL overlap
- SNPs: all overlapping nominal cis SNPs within the 1 Mb tensorQTL window
- minimum overlap: 10 SNPs per locus
- minimum eQTL evidence: at least one finite `abs(slope / slope_se) >= 2`
- GWAS filter: `SI >= 0.8`, no p-value cutoff
- method: `coloc.abf` plus `coloc::sensitivity`, not SuSiE
- parallelism: default `--n-cores 4`; set explicitly with `--n-cores N`

`--n-cores` sets the number of BiocParallel worker processes. It does not mean
OpenMP, BLAS, Arrow, or data.table threads per worker. `04_run_coloc.R` caps
those library threads to 1 before package loading so a run such as
`--n-cores 12` means 12 coloc workers, not 12 workers times many BLAS threads.

Per-disorder outputs are written under:

```text
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.qs2
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.runmeta.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.sensitivity.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.complete
```

Completed datasets are skipped only when the result, metadata, sensitivity, and
`.complete` marker all exist. If the marker is missing but the three primary
files exist, the script validates them and writes the marker before skipping. To
rerun a completed dataset, delete its existing coloc output files first.

## Summaries

`03_eqtl_explore.Rmd` defaults to:

```r
RUN_SEURAT <- TRUE
RUN_CUSTOM_CLUSTER <- FALSE
WRITE_DEG_VIEW_TABLES <- FALSE
```

The `03*.Rmd` notebooks in this folder process only the `all` group by default.
Nominal tables should also contain only `split == "all"` unless sex-stratified
analysis is explicitly re-enabled in a separate workflow.

Run the exact-tie recovery audit after `03_eqtl_explore.Rmd` and before using
significant eQTL variant IDs for boxplots or GWAS/coloc overlap interpretation:

```bash
Rscript ./03ab_recover_indepeQTLs.R --promote-current-pairs
```

The recovery starts from `map_significant_pairs.csv.gz`, uses the full nominal
parquet files to find variants with exactly equal nominal p-values, and then
requires donor-level genotype dosages to be identical up to allele orientation
(`abs(r) = 1` within numerical tolerance). This dosage check makes the recovery
valid for conditionally independent signals too: equivalent dosage vectors
remain equivalent after conditioning. The script preserves the reported signal
rank, expands each signal to all verified tied variants,
refreshes exact-variant GWAS annotations, and audits strong gated coloc lead-SNP
matches. With `--promote-current-pairs`, it preserves the original
single-variant-per-signal input separately and replaces
`map_significant_pairs.csv.gz` with the expanded table. The promoted canonical
table is the downstream source for exact variant-ID matching, the final eQTL
workbook, the final coloc lead-membership annotation, and future eQTL boxplot
selection. The parallel audit outputs retain the `tie_recovered` and
`tie_recovery` filename suffixes.
For an audit against staged snapshots, `--current-pairs=`, `--independent=`,
`--coloc-gated=`, and `--output-dir=` override those four paths without
modifying the primary files.

`05_coloc_explore.Rmd` loads coloc result objects and sensitivity tables,
writes flattened result, gated result, lead/SNP-level PP.H4, QC, and summary
tables, and plots strong coloc counts per domainCT by disorder.

Default outputs are written under:

- `processed-data/11_eQTL_coloc/seurat/tables`
- `processed-data/11_eQTL_coloc/seurat/coloc/tables`
- `plots/11_eQTL_coloc`

GWAS mixed-overlap details are documented in `README_GWAS-gene-lists.md`.

## DEG Support

The primary DEG overlap uses only the `broad_interaction` view.

Definition:

- Start from final DEG summary CSVs.
- Keep rows with F-test BH-adjusted `adj.P.Val < 0.05`.
- Take the union across PRECAST L-A, PRECAST L-R, Seurat L-A, and Seurat L-R.
- Match eQTL eGenes by `gene_id`.
- Do not context-match or sex-match.

Disorder-related DEG columns in `map_significant_summary.csv` use the same
F-test-filtered source rows, then require a same-disorder post-hoc t-test flag:

- MDD-related: `F_NTC.MDD_ttest` or `M_NTC.MDD_ttest` is significant.
- BD-related: `F_NTC.BD_ttest` or `M_NTC.BD_ttest` is significant.

These columns are disorder-related, not formal disorder-specific-only calls.
The broad DEG columns remain the four-file F-test union. See
`README_DEGs_by_disorder.md` for the statistical rationale.

DEG source files:

- `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

Use `gene_id` as the primary overlap key. Keep `gene_name` for reporting and validation.

Preferred wording:

- F-test significant dx/sex DGE genes
- sex-stratified diagnosis-associated DGE genes
- broad PRECAST-or-Seurat DGE-supported genes

Avoid claiming that every gene has a formal interaction-only effect.

## Environment

Python environment notes are in `README_venv.md`.

R/Rmd code in this folder should anchor paths with:

```r
here::i_am(".git/HEAD")
```
