# eQTL and Colocalization Workflow

This folder contains the MBv eQTL and colocalization workflow.

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
- `processed-data/ref/GWAS/{BPD,MDD,SCZD}`

Canonical Seurat eQTL paths:

- `processed-data/11_eQTL_coloc/seurat/tqtl_in`
- `processed-data/11_eQTL_coloc/seurat/tqtl_out`
- `processed-data/11_eQTL_coloc/seurat/tables`

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
- `utils.R`: shared DEG, GWAS, tensorQTL, summary, and plotting helpers.

## Prepare Inputs

Stage required local inputs:

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

Default outputs are written under:

- `processed-data/11_eQTL_coloc/seurat/tables`
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
- BPD-related: `F_NTC.BPD_ttest` or `M_NTC.BPD_ttest` is significant.

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
