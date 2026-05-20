# eQTL Derived Tables

Derived tables are written separately for each tensorQTL analysis series:

- `seurat/tables/`
- `custom_cluster/tables/`

Raw tensorQTL outputs remain in each series `tqtl_out/` directory.

## Default Tables

### `manifest_qc.csv`

Input preparation inventory copied from `tqtl_in/prep_manifest.csv`.

Columns:

- `dataset_id`: tensorQTL dataset prefix used for input/output files.
- `context`: annotation context represented by the dataset.
- `split`: donor subset: `all`, `male`, or `female`.
- `covariate_model`: covariate model used during tensorQTL input preparation.
- `n_samples`: number of donor pseudobulk samples.
- `n_genes_bed`: genes written to the tensorQTL expression BED.
- `n_genes_pca`: genes used for expression-PC estimation.
- `n_expr_pcs`: expression PCs added to the covariate file.

### `map_cis_summary.csv`

Context/split summary of significant lead cis-eQTL eGenes from tensorQTL `map_cis`, using `qval < 0.05`.

Columns:

- `n_eGenes`: number of significant eGenes.
- `n_SCZD_GWAS`: significant eGenes whose lead variant is SCZD GWAS-significant.
- `n_DEG`: significant eGenes in Jacqui's PRECAST+Seurat DEG union.
- `n_DEG_SCZD_GWAS`: significant eGenes that satisfy both DEG-union and SCZD GWAS support.
- `SCZD_GWAS_genes`: gene symbols for SCZD GWAS-supported eGenes.
- `DEG_genes`: gene symbols for DEG-union eGenes.
- `DEG_SCZD_GWAS_genes`: gene symbols satisfying both filters.

### `map_independent_summary.csv`

Context/split summary of significant independent cis-eQTL signals.

Independent signals are retained when the parent cis result has `qval < 0.05` and the independent signal has `pval_perm < 0.05`.

Columns:

- `n_independent_signals`: number of retained independent signals.
- `n_eGenes`: number of eGenes represented by those signals.
- `n_SCZD_GWAS`: retained independent signals whose variant is SCZD GWAS-significant.
- `n_DEG`: eGenes in Jacqui's PRECAST+Seurat DEG union.
- `n_DEG_SCZD_GWAS`: retained independent signals satisfying both DEG-union and SCZD GWAS support.
- `SCZD_GWAS_genes`: gene symbols with SCZD GWAS-supported independent signals.
- `DEG_genes`: DEG-union gene symbols represented by retained independent signals.
- `DEG_SCZD_GWAS_genes`: gene symbols satisfying both filters.

### `map_cis_significant.csv.gz`

Row-level significant lead cis-eQTL eGenes from tensorQTL `map_cis`, using `qval < 0.05`.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `DEG`: 1 if `gene_id` is in Jacqui's PRECAST+Seurat DEG union, otherwise 0.
- `SCZD_GWAS`: 1 if `variant_id` is SCZD GWAS-significant at the standard threshold, otherwise 0.
- `is_hla_gene`: TRUE when `gene_name` matches `^HLA-`.
- `is_mhc_gene`: TRUE when the gene overlaps hg38 `chr6:25000000-35000000`.
- `lead_variant_mhc`: TRUE when the lead variant is in hg38 `chr6:25000000-35000000`.
- `exclude_hla`: TRUE for HLA genes.
- `exclude_mhc`: TRUE when `is_mhc_gene` or `lead_variant_mhc` is TRUE.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`, `qval`, and `pval_nominal_threshold`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_cis_region_sensitivity_summary.csv`

Context/split summary of lead cis-eQTL eGenes under HLA/MHC sensitivity filters.

Columns:

- `region_filter`: `all`, `exclude_hla`, or `exclude_mhc`.
- Existing summary columns from `map_cis_summary.csv`.
- `n_HLA_eGenes`, `n_MHC_eGenes`, `n_HLA_DEG`, `n_MHC_DEG`, `n_lead_variant_MHC`.
- `HLA_genes`, `HLA_DEG_genes`, `MHC_genes`, `MHC_DEG_genes`.

### `map_independent_significant.csv.gz`

Row-level significant independent cis-eQTL signals.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `qval_parent`: parent `map_cis` q-value for the eGene.
- `DEG`: 1 if `gene_id` is in Jacqui's PRECAST+Seurat DEG union, otherwise 0.
- `SCZD_GWAS`: 1 if `variant_id` is SCZD GWAS-significant at the standard threshold, otherwise 0.
- `is_hla_gene`: TRUE when `gene_name` matches `^HLA-`.
- `is_mhc_gene`: TRUE when the gene overlaps hg38 `chr6:25000000-35000000`.
- `lead_variant_mhc`: TRUE when the independent-signal variant is in hg38 `chr6:25000000-35000000`.
- `exclude_hla`: TRUE for HLA genes.
- `exclude_mhc`: TRUE when `is_mhc_gene` or `lead_variant_mhc` is TRUE.

Rows are retained when `qval_parent < 0.05` and `pval_perm < 0.05`.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`, and `rank`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_independent_region_sensitivity_summary.csv`

Context/split summary of independent cis-eQTL signals under HLA/MHC sensitivity filters.

Columns:

- `region_filter`: `all`, `exclude_hla`, or `exclude_mhc`.
- Existing summary columns from `map_independent_summary.csv`.
- `n_HLA_eGenes`, `n_MHC_eGenes`, `n_HLA_DEG`, `n_MHC_DEG`, `n_lead_variant_MHC`.
- `HLA_genes`, `HLA_DEG_genes`, `MHC_genes`, `MHC_DEG_genes`.

HLA/MHC sensitivity summaries are intended for interpretation and reporting. HLA/MHC colocalization results should be treated as exploratory unless a dedicated HLA-aware analysis is performed, because this region has unusual polymorphism, gene density, and linkage disequilibrium.

## Optional DEG-View Tables

The summary code has a disabled-by-default guard for DEG-view tables. These context/sex-localized breakdowns are not part of the default table set.
