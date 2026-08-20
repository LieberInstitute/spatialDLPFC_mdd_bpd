# Extra Models and DEG Views

This file keeps secondary context for eQTL models that are not the default workflow in `README.md`.

Default workflow remains Seurat clusters and all donors. Its `DEG` annotation is
the authoritative 503-gene set from `code/09_DEG_GRN/load_DEGs.r::sig.df`.

## Project Context

The MBv project studies sex-specific gene expression changes associated with MDD and BPD across dlPFC cortical layers.

Data summary:

- 119 adult postmortem human brain donors.
- Diagnostic groups: NTC, MDD, and BPD.
- 10x Genomics Visium spatial RNA-seq.
- PRECAST/smoothed annotation: 6 spatial domains, `L1`, `L2`, `L3.4`, `L5`, `L6`, `WM`.
- Seurat label-transfer annotation: 8 cell-type-like clusters, `Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`, `Oligo`.

The eQTL workflow uses donor-level pseudobulk expression and genotype data. The current primary model avoids genotype-by-diagnosis interaction eQTLs and instead maximizes donor count for cis-eQTL discovery.

## DGE Model Context

The DGE analyses use `limma-voom` models on donor-level pseudobulk expression with diagnosis, sex, annotation context, global PC3, age, and `nspots` covariates.

The final `dx-sex_degs-F-test-t-test.csv` files are generated from sex-stratified diagnosis contrasts such as `F_NTC.MDD`, `M_NTC.MDD`, `F_NTC.BPD`, and `M_NTC.BPD`, followed by an omnibus F-test across tested contrasts.

Recommended wording:

- F-test significant dx/sex DGE genes
- sex-stratified diagnosis-associated DGE genes
- broad dx/sex DGE-supported genes

Avoid: all genes with a formal diagnosis-by-sex interaction-only effect.

Model families:

1. Layer-/context-adjusted (L-A): `~ 0 + condition:sex + context + PC3 + age + nspots`
2. Layer-/context-restricted (L-R): `~ 0 + condition:context:sex + PC3 + age + nspots`

In file names, `layer` should be read as annotation context:

- PRECAST/smoothed runs use spatial domains.
- Seurat runs use cell-type-like clusters.

## Authoritative DEG Universe

The primary project DEG universe is `unique(sig.df$gene_id)` after sourcing
`code/09_DEG_GRN/load_DEGs.r`. The loader combines the following four omnibus
F-test summaries with their matching `*_rev-gene-input_moderated-t-test.csv`
files and retains `adj.P.Val2 < 0.05`:

- PRECAST L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- PRECAST L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- Seurat L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- Seurat L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

This produces exactly 503 one-to-one gene ID/name pairs. The earlier direct
union of the four F-test summaries contains 817 genes and remains available only
through `load_legacy_standard_DEGs()` for historical comparison.

Historical project notes describe the legacy broad union as useful because the
Seurat eQTL pseudobulks are not identical to either PRECAST DGE inputs or Seurat
DGE inputs:

- PRECAST excludes low-UMI spots.
- Seurat DEG models include low-UMI cluster spots.
- Current Seurat eQTL pseudobulks aggregate by cell-type labels while excluding low-UMI spots.

## DEG Views

The overlap code can report four DEG views derived from authoritative `sig.df`.
Normal final-table annotation uses the global 503-gene set and does not
context-match or sex-match it.

| DEG view | DEG files used | Extra DEG filter after `adj.P.Val < 0.05` | eQTL matching |
| --- | --- | --- | --- |
| `broad_interaction` | all authoritative `sig.df` rows | none | gene only |
| `context_localized` | Seurat `L-A` plus matching Seurat context rows | significant `sig.df` rows | gene and Seurat context |
| `sex_specific` | all authoritative `sig.df` rows | matching `F_` or `M_` contrast | gene and sex |
| `context_and_sex_specific` | matching Seurat context rows | matching context and `F_` or `M_` contrast | gene, Seurat context, and sex |

All views start from significant `sig.df` rows, which already combine omnibus
F-test eligibility with adjusted moderated t-test significance.

### broad_interaction

Gene-level set of the 503 significant genes in `sig.df`.

Matching:

- not context-matched
- not sex-matched
- gene-level overlap only

### context_localized

Source:

- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

Filter:

- same Seurat context has `n_ttest_sig_<context> > 0`

Matching:

- context-matched
- not sex-matched

This view excludes PRECAST-only DGE support because PRECAST domains do not map one-to-one onto Seurat eQTL contexts.

### sex_specific

Source:

- PRECAST L-A
- PRECAST L-R
- Seurat L-A
- Seurat L-R

Filter:

- female eQTLs use any `F_*_ttest`
- male eQTLs use any `M_*_ttest`

Matching:

- sex-matched
- not context-matched

### context_and_sex_specific

Source:

- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

Filter:

- same context and same sex-specific t-test evidence
- example: female Astro eQTLs use columns such as `Astro_F_NTC.MDD_ttest`

Matching:

- context-matched
- sex-matched

This is the strictest DEG overlap view and is expected to have smaller counts.

## Sex-Stratified eQTLs

Sex-stratified eQTL inputs are optional.

Prepare all Seurat splits:

```bash
Rscript ./01_prep_inputs.R --splits all,male,female
```

Prepare one sex split:

```bash
Rscript ./01_prep_inputs.R --splits female
Rscript ./01_prep_inputs.R --splits male
```

Run all Seurat splits:

```bash
./02_run_tensorQTL.sh --include-sex-splits
```

Run selected splits:

```bash
./02_run_tensorQTL.sh --splits female
./02_run_tensorQTL.sh --splits male
```

Dataset suffixes:

- all donors: no suffix, for example `astro`
- male: `_m`, for example `astro_m`
- female: `_f`, for example `astro_f`

Covariates:

- all donors: `DX + sex + age + PC3 + nspots + snpPC1-5 + expression PCs`
- male/female: `DX + age + PC3 + nspots + snpPC1-5 + expression PCs`

To include sex splits in `03_eqtl_explore.Rmd`, set:

```r
INCLUDE_SEX_SPLITS <- TRUE
```

## Custom-Cluster Inputs

The custom-cluster workflow is optional and does not overwrite Seurat inputs.

Custom input object:

- `processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata`

Primary custom model:

- grouping column: `custom_cluster`
- split: `all`
- covariates: `DX + sex + age + PC3 + nspots + snpPC1-5 + expression PCs`
- input directory: `processed-data/11_eQTL_coloc/custom_cluster/tqtl_in`
- output directory: `processed-data/11_eQTL_coloc/custom_cluster/tqtl_out`

Stage custom-cluster inputs:

```bash
./stage_required_data.sh --include-custom-cluster
```

Prepare custom-cluster tensorQTL inputs:

```bash
Rscript ./01_prep_inputs.R \
  --cluster-col custom_cluster \
  --analysis-label custom_cluster \
  --model nspots \
  --splits all \
  --out-dir ../../processed-data/11_eQTL_coloc/custom_cluster/tqtl_in
```

Run custom-cluster tensorQTL:

```bash
./02_run_tensorQTL.sh --analysis custom_cluster
```

Dry-run custom-cluster commands:

```bash
./02_run_tensorQTL.sh --analysis custom_cluster --dry-run
```

To include custom clusters in `03_eqtl_explore.Rmd`, set:

```r
RUN_CUSTOM_CLUSTER <- TRUE
```

## Custom-Cluster Construction Notes

The custom-cluster strategy combines Seurat and PRECAST labels. It starts from Seurat labels, removes low-UMI spots, and uses PRECAST domains to split ambiguous labels such as Astro and L2.3.

Example custom labels:

- `Astro.L1`
- `Astro.Nrn`
- `Inhb`
- `L2`
- `L3`
- `L4`
- `L5`
- `L6`
- `Micro.Vasc`
- `WM`

Related source paths:

- `code/06_pseudobulk/custom_cluster/01_create_pseudobulk.r`
- `code/06_pseudobulk/custom_cluster/02_norm-QC-PCA.r`

## Optional Rmd Controls

`03_eqtl_explore.Rmd` has these top-level controls:

```r
RUN_SEURAT <- TRUE
RUN_CUSTOM_CLUSTER <- FALSE
INCLUDE_SEX_SPLITS <- FALSE
WRITE_DEG_VIEW_TABLES <- FALSE
```

Set `WRITE_DEG_VIEW_TABLES <- TRUE` to emit DEG-view summary tables.
