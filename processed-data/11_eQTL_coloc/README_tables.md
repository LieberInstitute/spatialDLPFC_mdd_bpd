# eQTL Derived Tables

Derived tables are written separately for each tensorQTL analysis series:

- `seurat/tables/`
- `custom_cluster/tables/`

Raw tensorQTL outputs remain in each series `tqtl_out/` directory when they are
staged locally. The current shared `seurat/tables/` snapshot is usable for
review, but it has a mixed column schema: the summary tables and
`eQTL_boxplot_deg_pairs.csv` use the cleaned GWAS column names, while several
row-level `.csv.gz` files still use the older internal `*_GWAS` / `*_GWASg`
names. This README documents both so table consumers can interpret the current
files without regenerating them.

## Which Tables To Use

- Use `map_significant_summary.csv` for downstream integrated eGene counts by
  context and GWAS set.
- Use `map_significant_pairs.csv.gz` when one row per significant
  gene-variant pair is needed.
- Use `map_significant_unified.csv.gz` when source-specific cis versus
  independent tensorQTL statistics must be audited.
- Use `eQTL_boxplot_deg_pairs.csv` only as genotype-boxplot metadata for
  broad-DEG significant eQTL pairs. It is not an all-eGene or all-pair summary.
- Treat `map_cis_summary.csv` and `map_independent_summary.csv` as summaries of
  the intermediate cis and independent result sets used to build the final
  significant-pair union.

## GWAS Overlap Conventions

Disorder keys in cleaned output columns are `MDD`, `BD`, and `SCZD`. Older
row-level tables may still use `BPD`; that refers to the same bipolar disorder
GWAS set now labeled `BD`.

Strict GWAS overlap uses genome-wide significant variants. Exploratory overlap
is a mood-disorder-only exact-variant check for MDD and BD at the relaxed
suggestive threshold used in the current workflow.

Cleaned row-level columns:

- `DIS_gwasVar_strict`: `1` when the row's eQTL `variant_id` is in the strict
  GWAS variant set for `DIS`; otherwise `0`.
- `DIS_gwasVar_exp`: `1` when the row's eQTL `variant_id` is in the
  exploratory GWAS variant set for `DIS`; otherwise `0`. Present for MDD and BD
  only.
- `DIS_gwasGene`: `1` when the row's `gene_name` is in the curated GWAS gene
  list for `DIS`; otherwise `0`.
- `DIS_gwas_strict`: combined strict support, equal to
  `DIS_gwasVar_strict OR DIS_gwasGene`.
- `DIS_gwas_exp`: combined exploratory support, equal to
  `DIS_gwasVar_exp OR DIS_gwasGene`. Present for MDD and BD only.
- `DIS_gwasP`, `DIS_gwasBeta`, `DIS_gwasBetaSE`: exploratory exact-variant
  GWAS statistics for MDD/BD rows with `DIS_gwasVar_exp == 1`; otherwise
  missing.

Cleaned summary columns count unique eGenes within each `context`/`split` row:

- `n_DIS_gwas_strict`: unique eGenes with combined strict support
  (`gwasVar_strict OR gwasGene`), not a variant-only count.
- `n_DIS_gwas_exp`: unique eGenes with combined exploratory support
  (`gwasVar_exp OR gwasGene`), not a variant-only count.
- `DIS_gwas_strict_genes` and `DIS_gwas_exp_genes`: gene symbols contributing
  to those combined counts.
- `n_trifecta_DIS_strict` / `n_trifecta_DIS_exp`: unique broad-DEG eGenes with
  combined eQTL + DEG + GWAS support.
- `n_DIS_DEGxGWAS_strict` / `n_DIS_DEGxGWAS_exp`: unique disorder-related DEG
  eGenes with combined same-disorder GWAS support.

The current `map_significant_summary.csv` does not expose dedicated variant-only
summary columns such as `n_DIS_gwasVar_strict`. For exact variant-overlap counts,
use a row-level table that contains `DIS_gwasVar_*` columns or regenerate a
summary after the required raw inputs are staged.

Older/internal row-level columns:

- `DIS_GWAS`: strict exact-variant support.
- `DIS_GWASg_gene`: curated GWAS gene-list support.
- `DIS_GWASg_variant`: duplicate of `DIS_GWAS`.
- `DIS_GWASg`: combined strict support, equal to
  `DIS_GWAS OR DIS_GWASg_gene`.

Older/internal row-level files in the current snapshot do not contain the
exploratory `*_exp` columns.

## Default Tables

### `manifest_qc.csv`

Input preparation inventory copied from `tqtl_in/prep_manifest.csv`.

Columns:

- `dataset_id`: tensorQTL dataset prefix used for input/output files.
- `context`: annotation context represented by the dataset.
- `split`: donor subset, usually `all` for the current default Seurat analysis.
- `covariate_model`: covariate model used during tensorQTL input preparation.
- `n_samples`: number of donor pseudobulk samples.
- `n_genes_bed`: genes written to the tensorQTL expression BED.
- `n_genes_pca`: genes used for expression-PC estimation.
- `n_expr_pcs`: expression PCs added to the covariate file.

### `map_cis_summary.csv`

Context/split summary of the intermediate significant lead cis-eQTL eGene set
from tensorQTL `map_cis`, using `qval < 0.05`.

This table is useful for auditing the lead-cis input to the final union, but the
preferred downstream integrated eGene summary is `map_significant_summary.csv`.

Common columns:

- `split`, `context`: grouping variables.
- `n_eGenes`: number of significant lead cis-eQTL eGenes.
- `n_DEG`: number of significant lead cis-eQTL eGenes in the broad
  PRECAST-or-Seurat DEG union.
- `n_DIS_gwas_strict`: unique significant lead cis-eQTL eGenes with combined
  strict GWAS support for `DIS`.
- `n_trifecta_DIS_strict`: unique broad-DEG lead cis-eQTL eGenes with combined
  strict GWAS support for `DIS`.
- `DEG_genes`, `DIS_gwas_strict_genes`, `trifecta_DIS_strict_genes`: gene
  symbol lists for the corresponding counts.

### `map_independent_summary.csv`

Context/split summary of the intermediate significant independent cis-eQTL
signal set.

Independent signals are retained when the parent cis result has `qval < 0.05`
and the independent signal has `pval_perm < 0.05`. This table is useful for
auditing the independent-signal input to the final union; the preferred
downstream integrated eGene summary is `map_significant_summary.csv`.

Common columns:

- `n_independent_signals`: number of retained independent signals.
- `n_eGenes`: number of eGenes represented by those signals.
- `n_DEG`: number of represented eGenes in the broad PRECAST-or-Seurat DEG
  union.
- `n_DIS_gwas_strict`: unique independent-signal eGenes with combined strict
  GWAS support for `DIS`.
- `n_trifecta_DIS_strict`: unique broad-DEG independent-signal eGenes with
  combined strict GWAS support for `DIS`.
- `DEG_genes`, `DIS_gwas_strict_genes`, `trifecta_DIS_strict_genes`: gene
  symbol lists for the corresponding counts.

### `map_significant_summary.csv`

Context/split summary of `map_significant_pairs.csv.gz`, the simplified
one-row-per-significant-pair eQTL table. This is the recommended integrated
summary for eGene-centered downstream reporting.

Rows are context/split summaries, not row-level eQTL pairs.

Columns:

- `n_significant_pairs`: number of retained unique gene-variant pairs.
- `n_cis_supported_pairs`: number of retained pairs present in
  `map_cis_significant.csv.gz`.
- `n_indep_supported_pairs`: number of retained pairs present in
  `map_independent_significant.csv.gz`.
- `n_shared_pairs`: number of retained pairs present in both source tables.
- `n_cis_only_pairs`: number of retained pairs present only in
  `map_cis_significant.csv.gz`.
- `n_indep_only_pairs`: number of retained pairs present only in
  `map_independent_significant.csv.gz`.
- `n_eGenes`: number of unique eGenes represented by retained pairs.
- `n_DEG`: number of retained eGenes in the broad PRECAST-or-Seurat DEG union.
- `n_MDD_DEG`, `n_BD_DEG`: retained eGenes in the same-disorder
  contrast-supported DEG sets.
- `n_DIS_gwas_strict`, `n_DIS_gwas_exp`: combined GWAS-supported eGene counts as
  described in "GWAS Overlap Conventions".
- `n_trifecta_DIS_strict`, `n_trifecta_DIS_exp`: broad-DEG + eQTL + combined
  GWAS-supported eGene counts.
- `n_DIS_DEGxGWAS_strict`, `n_DIS_DEGxGWAS_exp`: same-disorder DEG + eQTL +
  combined same-disorder GWAS-supported eGene counts.
- `*_genes` columns: comma-separated gene-symbol lists matching the count
  columns with the same prefix.

`n_significant_pairs` can exceed `n_eGenes` when an eGene has multiple retained
variants, including conditionally independent secondary signals.

### `map_cis_significant.csv.gz`

Row-level significant lead cis-eQTL eGenes from tensorQTL `map_cis`, using
`qval < 0.05`.

The current shared snapshot uses the older/internal strict-only GWAS columns in
this file. See "GWAS Overlap Conventions" for the mapping from
`DIS_GWAS`/`DIS_GWASg` to cleaned meanings.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `gene_id`, `gene_name`, `variant_id`
- `DEG`: `1` if `gene_id` is in the broad PRECAST-or-Seurat DEG union,
  otherwise `0`.
- Older/internal strict GWAS support columns such as `MDD_GWAS`, `BPD_GWAS`,
  `SCZD_GWAS`, and `*_GWASg`.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`,
`qval`, and `pval_nominal_threshold`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_independent_significant.csv.gz`

Row-level significant independent cis-eQTL signals.

Rows are retained when `qval_parent < 0.05` and `pval_perm < 0.05`. The current
shared snapshot uses the older/internal strict-only GWAS columns in this file.

The table preserves tensorQTL-native columns and adds project annotations:

- `dataset_id`, `context`, `split`
- `gene_id`, `gene_name`, `variant_id`
- `qval_parent`: parent `map_cis` q-value for the eGene.
- `DEG`: `1` if `gene_id` is in the broad PRECAST-or-Seurat DEG union,
  otherwise `0`.
- Older/internal strict GWAS support columns such as `MDD_GWAS`, `BPD_GWAS`,
  `SCZD_GWAS`, and `*_GWASg`.

For tensorQTL-native columns such as `num_var`, `beta_shape1`, `pval_perm`, and
`rank`, see:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

### `map_significant_unified.csv.gz`

Lossless long-format union of `map_cis_significant.csv.gz` and
`map_independent_significant.csv.gz`.

Rows are source-specific observations, not unique gene-variant pairs. This
preserves all tensorQTL-native values because the same gene-variant pair can
have different source-specific values for columns such as `pval_nominal`,
`slope`, `slope_se`, `pval_perm`, `pval_beta`, and beta-shape fields.

Added columns:

- `result_source`: `cis` for rows from `map_cis_significant.csv.gz`, or
  `independent` for rows from `map_independent_significant.csv.gz`.
- `pair_provenance`: `cis_only`, `independent_only`, or `cis_and_independent`,
  computed by `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, and
  `variant_id`.

The current shared snapshot uses older/internal strict-only GWAS columns in this
file.

### `map_significant_pairs.csv.gz`

Simplified one-row-per-pair table derived from `map_significant_unified.csv.gz`.

Rows are unique by `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, and
`variant_id`. The `source` column identifies which source row supplied the
retained tensorQTL statistics:

- `cis`: the pair was present in `map_cis_significant.csv.gz`; this is also
  used for shared cis+independent pairs.
- `indep`: the pair was present only in `map_independent_significant.csv.gz`.

Support/provenance columns:

- `pair_provenance`: `cis_only`, `independent_only`, or `cis_and_independent`.
- `cis_supported`: `TRUE` if the pair was present in `map_cis_significant.csv.gz`.
- `indep_supported`: `TRUE` if the pair was present in
  `map_independent_significant.csv.gz`.

The table keeps the selected source row's tensorQTL statistic column names
unchanged. It intentionally does not preserve duplicate source-specific
statistics for shared pairs; use `map_significant_unified.csv.gz` for that audit
trail.

The current shared snapshot uses older/internal strict-only GWAS columns in this
file. Its strict combined counts match the strict combined counts exposed in
`map_significant_summary.csv`.

### `eQTL_boxplot_deg_pairs.csv`

Broad-DEG significant eQTL pair metadata used by `03b_eQTL_boxplots.Rmd` to make
genotype boxplots.

Rows are selected from significant eQTL pairs with `split == "all"` and
`DEG == 1`. The current file has 263 rows and is unique by
`dataset_id + gene_id + variant_id`. It contains one retained plotting row per
broad-DEG significant pair after ordering by context, gene, nominal p-value, and
variant.

Use this table for boxplot pair metadata, not for all-eGene or all-significant
eQTL counts. It is a filtered plotting table, so its counts are expected to be
smaller than `map_significant_summary.csv`.

Important columns:

- `example_type`, `example_deg_basis`: provenance for the plotting subset.
  Current rows are `broad_deg_significant` / `broad_DEG`.
- `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, `phenotype_id`,
  `variant_id`: dataset, context, and eQTL pair identifiers.
- `DEG`: broad PRECAST-or-Seurat DEG union flag; all current rows have
  `DEG == 1`.
- `MDD_DEG`, `BD_DEG`: same-disorder contrast-supported DEG flags.
- `result_source`, `pair_provenance`, `cis_supported`, `indep_supported` when
  present: significant-pair provenance inherited from the final eQTL union.
- `DIS_gwasVar_strict`, `DIS_gwasVar_exp`, `DIS_gwasGene`,
  `DIS_gwas_strict`, `DIS_gwas_exp`, `DIS_gwasP`, `DIS_gwasBeta`,
  `DIS_gwasBetaSE`: cleaned GWAS overlap columns described in
  "GWAS Overlap Conventions".
- `rsid`, `rsid_source`: dbSNP label and source from the genotype-wide variant
  cache.

The combined flags in this file are internally consistent:
`DIS_gwas_strict == DIS_gwasVar_strict OR DIS_gwasGene`, and for MDD/BD
`DIS_gwas_exp == DIS_gwasVar_exp OR DIS_gwasGene`.

### `nominal_BH05.csv.gz` and `nominal_BH05.xlsx`

BH-adjusted significant nominal eQTL pairs from tensorQTL `map_nominal`, using
`fdr < 0.05` within each dataset/context. The `.xlsx` workbook mirrors the
nominal table plus summary and QC sheets for review.

Columns:

- `dataset_id`, `context`, `split`, `gene_id`, `gene_name`, `phenotype_id`,
  `variant_id`: dataset and pair identifiers.
- `pval_nominal`, `fdr`, `slope`, `slope_se`, `start_distance`, `af`,
  `ma_samples`, `ma_count`: nominal eQTL statistics from tensorQTL.
- `DEG` and GWAS overlap columns: project annotations using either the cleaned
  or older/internal schema depending on the specific file snapshot.

### `nominal_BH05_summary.csv`

Context/split summary of `nominal_BH05.csv.gz`.

Columns:

- `n_nominal_pairs`: number of retained BH-significant nominal gene-variant
  pairs.
- `n_eGenes`: number of eGenes represented by those pairs.
- DEG and GWAS count/gene-list columns use the same unique-eGene definitions as
  the other summary tables.

### `nominal_BH05_qc.csv`

QC inventory for nominal eQTL aggregation.

Columns:

- `n_parquet_files`: per-chromosome tensorQTL nominal parquet files read.
- `n_nominal_rows_tested`: nominal pairs tested before filtering.
- `n_p001_rows`: nominal pairs with `pval_nominal <= 0.001`.
- `n_BH05_rows`: nominal pairs retained at `fdr < 0.05`.
- `n_BH05_genes`: genes represented by retained nominal pairs.

## Regeneration Notes

Do not regenerate tables just to resolve documentation ambiguity. A full clean
regeneration of the current Seurat table set requires staged raw inputs that are
not present in the current local snapshot, including:

- `processed-data/11_eQTL_coloc/seurat/tqtl_in/prep_manifest.csv`
- `processed-data/11_eQTL_coloc/seurat/tqtl_out/*.gene.map_cis.tab.gz`
- `processed-data/11_eQTL_coloc/seurat/tqtl_out/*.gene.map_independent.tsv.gz`
- `processed-data/00_genotypes/plink2/merged_maf05.pvar`
- current `BD`, `MDD`, and `SCZD` GWAS caches or BCFs matching the code's
  strict and exploratory paths.

Only regenerate row-level `map_significant_*` or `eQTL_boxplot_deg_pairs.csv`
after those inputs are staged and a targeted audit shows that a table value is
wrong, not merely that the column names are confusing.

## Optional DEG-View Tables

The summary code has a disabled-by-default guard for DEG-view tables. These
context/sex-localized breakdowns are not part of the default table set.
