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

Current mood-disorder GWAS inputs include the supplied European 23andMe
component:

- BD: `processed-data/ref/GWAS/BD/bip2024_eur.hg38.bcf`, release tag
  `full23andMe_preDENTIST`.
- MDD: `processed-data/ref/GWAS/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf`,
  release tag `full23andMe`.

The integrated BD and MDD BCFs have no combined imputation-quality statistic.
Their cache names therefore use `SInone`; `SI >= 0.8` remains active for SCZD.
The reconstructed BD file remains pre-DENTIST and is not the exact final paper
summary-statistics file.

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
- `03a_nominal_eQTLs.Rmd`: gather nominal all-donor Seurat eQTLs.
- `03b_eQTL_boxplots.Rmd`: select example eQTL pairs and render genotype boxplots.
- `03c_GWAS_relaxed_eQTLs.R`: rebuild exploratory MDD/BD `p < 1e-5` exact-variant overlap tables.
- `04_run_coloc.R`: run coloc ABF and sensitivity checks from full tensorQTL nominal parquet.
- `05_coloc_explore.Rmd`: flatten coloc outputs, apply sensitivity gates, write tables, and plot strong coloc counts.
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
- GWAS filter: no p-value cutoff; `SI >= 0.8` where SI is available, and no SI
  filter for integrated BD/MDD because combined SI is unavailable
- method: `coloc.abf` plus `coloc::sensitivity`, not SuSiE
- parallelism: default `--n-cores 4`; set explicitly with `--n-cores N`

`--n-cores` sets the number of BiocParallel worker processes. It does not mean
OpenMP, BLAS, Arrow, or data.table threads per worker. `04_run_coloc.R` caps
those library threads to 1 before package loading so a run such as
`--n-cores 12` means 12 coloc workers, not 12 workers times many BLAS threads.

Per-disorder outputs are written under:

```text
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.runmeta.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.sensitivity.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.flat.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.errors.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/<DIS>/coloc_<dataset_id>.complete
```

Full `coloc_<dataset_id>.qs2` objects are written only with `--write-full-qs2`.

### 2026-07-28 integrated-GWAS rerun

The `new-gwas` branch contains a complete local rerun of `coloc::coloc.abf`
and `coloc::sensitivity` for all eight BD and MDD cell contexts using the
23andMe-inclusive GWAS inputs above. SCZD was not rerun. The new per-context
outputs and release-specific dense GWAS caches are under:

```text
processed-data/11_eQTL_coloc/seurat/coloc/{BD,MDD}/
```

The rerun used the lean default output mode, so it wrote candidate summaries,
sensitivity results, error tables, run metadata, and completion markers, but
not full `.qs2` result objects. Aggregate and final outputs are under:

```text
processed-data/11_eQTL_coloc/seurat/coloc/tables/
processed-data/11_eQTL_coloc/seurat/final/
```

Execution logs, the frozen no-23andMe results, and machine-readable comparison
tables are under:

```text
processed-data/11_eQTL_coloc/seurat/run_logs/23andMe_2026-07-28/
processed-data/11_eQTL_coloc/seurat/archive/no23andMe_2026-07-28/
processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28/
```

The comparison and GWAS provenance are documented in
`GWAS-23andMe-MDD-BD-impact-on-tables.md`.

Generated files retain their canonical names on both `devel` and `new-gwas`.
The branch is the version identifier: `devel` contains the prior tracked
workbooks and figures, while `new-gwas` contains their regenerated versions at
the same paths. Temporary `*-new-gwas.<ext>` copies were moved back to those
canonical paths before commit and therefore no longer exist. Use the dated
archive and comparison directories for side-by-side local files rather than
adding filename suffixes to workflow outputs.

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

`05_coloc_explore.Rmd` loads coloc result objects and sensitivity tables,
writes flattened result, gated result, lead/SNP-level PP.H4, QC, and summary
tables, and plots strong coloc counts per domainCT by disorder.

Default outputs are written under:

- `processed-data/11_eQTL_coloc/seurat/tables`
- `processed-data/11_eQTL_coloc/seurat/coloc/tables`
- `plots/11_eQTL_coloc`

GWAS mixed-overlap details are documented in `README_GWAS_exploratory.md` and
`processed-data/ref/GWAS/README_genes_variants_lists.md`.

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
