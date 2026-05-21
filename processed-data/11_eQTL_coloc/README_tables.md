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
- `n_<DIS>_GWAS`: significant eGenes whose lead variant is GWAS-significant for `<DIS>`.
- `n_<DIS>_GWASg`: significant eGenes with either GWAS-significant variant support or exact GWAS gene-list support for `<DIS>`.
- `n_DEG`: significant eGenes in Jacqui's PRECAST+Seurat DEG union.
- `n_DEG_<DIS>_GWAS`: significant eGenes that satisfy both DEG-union and variant-only GWAS support.
- `n_DEG_<DIS>_GWASg`: significant eGenes that satisfy both DEG-union and mixed GWASg support. Use this for primary DEG x eQTL x GWAS reporting.
- `<DIS>_GWAS_genes`: gene symbols for variant-only GWAS-supported eGenes.
- `<DIS>_GWASg_genes`: gene symbols for mixed GWASg-supported eGenes.
- `DEG_genes`: gene symbols for DEG-union eGenes.
- `DEG_<DIS>_GWAS_genes`: gene symbols satisfying both DEG and variant-only GWAS filters.
- `DEG_<DIS>_GWASg_genes`: gene symbols satisfying both DEG and mixed GWASg filters.

### `map_independent_summary.csv`

Context/split summary of significant independent cis-eQTL signals.

Independent signals are retained when the parent cis result has `qval < 0.05` and the independent signal has `pval_perm < 0.05`.

Columns:

- `n_independent_signals`: number of retained independent signals.
- `n_eGenes`: number of eGenes represented by those signals.
- `n_<DIS>_GWAS`: retained independent-signal eGenes whose variant is GWAS-significant for `<DIS>`.
- `n_<DIS>_GWASg`: retained independent-signal eGenes with either GWAS-significant variant support or exact GWAS gene-list support for `<DIS>`.
- `n_DEG`: eGenes in Jacqui's PRECAST+Seurat DEG union.
- `n_DEG_<DIS>_GWAS`: retained independent-signal eGenes satisfying both DEG-union and variant-only GWAS support.
- `n_DEG_<DIS>_GWASg`: retained independent-signal eGenes satisfying both DEG-union and mixed GWASg support. Use this for primary DEG x eQTL x GWAS reporting.
- `<DIS>_GWAS_genes`: gene symbols with variant-only GWAS-supported independent signals.
- `<DIS>_GWASg_genes`: gene symbols with mixed GWASg-supported independent signals.
- `DEG_genes`: DEG-union gene symbols represented by retained independent signals.
- `DEG_<DIS>_GWAS_genes`: gene symbols satisfying both DEG and variant-only GWAS filters.
- `DEG_<DIS>_GWASg_genes`: gene symbols satisfying both DEG and mixed GWASg filters.

### `map_significant_summary.csv`

Context/split summary of `map_significant_pairs.csv.gz`, the simplified one-row-per-pair eQTL table.

Columns:

- `n_significant_pairs`: number of retained unique gene-variant pairs.
- `n_cis_pairs`: number of retained pairs whose statistics come only from `map_cis`.
- `n_indep_pairs`: number of retained pairs whose statistics come from `map_independent`.
- `n_eGenes`: number of eGenes represented by retained pairs.
- DEG, GWAS, and GWASg count/gene-list columns use the same unique-eGene definitions as the cis and independent summaries.

`n_eGenes` is expected to match the cis and independent summaries. `n_significant_pairs` can exceed `n_eGenes` when secondary independent signals exist.

### `map_cis_significant.csv.gz`

Row-level significant lead cis-eQTL eGenes from tensorQTL `map_cis`, using `qval < 0.05`.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `DEG`: 1 if `gene_id` is in Jacqui's PRECAST+Seurat DEG union, otherwise 0.
- `<DIS>_GWAS`: 1 if `variant_id` is GWAS-significant for `<DIS>` at the standard threshold, otherwise 0.
- `<DIS>_GWASg`: 1 if either `<DIS>_GWASg_variant` or `<DIS>_GWASg_gene` is 1.
- `<DIS>_GWASg_variant`: same variant-only evidence as `<DIS>_GWAS`.
- `<DIS>_GWASg_gene`: 1 if `gene_name` exactly matches the selected broad GWAS gene list for `<DIS>`.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`, `qval`, and `pval_nominal_threshold`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_independent_significant.csv.gz`

Row-level significant independent cis-eQTL signals.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `qval_parent`: parent `map_cis` q-value for the eGene.
- `DEG`: 1 if `gene_id` is in Jacqui's PRECAST+Seurat DEG union, otherwise 0.
- `<DIS>_GWAS`: 1 if `variant_id` is GWAS-significant for `<DIS>` at the standard threshold, otherwise 0.
- `<DIS>_GWASg`: 1 if either `<DIS>_GWASg_variant` or `<DIS>_GWASg_gene` is 1.
- `<DIS>_GWASg_variant`: same variant-only evidence as `<DIS>_GWAS`.
- `<DIS>_GWASg_gene`: 1 if `gene_name` exactly matches the selected broad GWAS gene list for `<DIS>`.

Rows are retained when `qval_parent < 0.05` and `pval_perm < 0.05`.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`, and `rank`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_significant_unified.csv.gz`

Lossless long-format union of `map_cis_significant.csv.gz` and `map_independent_significant.csv.gz`.

Rows are source-specific observations, not unique gene-variant pairs. This preserves all tensorQTL-native values because the same gene-variant pair can have different source-specific values for columns such as `pval_nominal`, `slope`, `slope_se`, `pval_perm`, `pval_beta`, and beta-shape fields.

Added columns:

- `result_source`: `cis` for rows from `map_cis_significant.csv.gz`, or `independent` for rows from `map_independent_significant.csv.gz`.
- `pair_provenance`: `cis_only`, `independent_only`, or `cis_and_independent`, computed by `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, and `variant_id`.

`map_cis` and `map_independent` overlap at eGene level, but not always at gene-variant-pair level. `map_cis` reports the lead cis-eQTL per phenotype from the permutation/q-value path. `map_independent` reports tensorQTL stepwise independent signals; rank 1 is usually, but not always, the same variant as `map_cis`, and ranks 2/3 can add conditionally independent secondary signals.

### `map_significant_pairs.csv.gz`

Simplified one-row-per-pair table derived from `map_significant_unified.csv.gz`.

Rows are unique by `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, and `variant_id`. The `source` column identifies which source row supplied the retained tensorQTL statistics:

- `indep`: the pair was present in `map_independent_significant.csv.gz`; this is also used for shared cis+independent pairs.
- `cis`: the pair was present only in `map_cis_significant.csv.gz`.

The table keeps the selected source row's tensorQTL statistic column names unchanged. It intentionally does not preserve duplicate source-specific statistics for shared pairs; use `map_significant_unified.csv.gz` for that audit trail.

## Optional DEG-View Tables

The summary code has a disabled-by-default guard for DEG-view tables. These context/sex-localized breakdowns are not part of the default table set.
