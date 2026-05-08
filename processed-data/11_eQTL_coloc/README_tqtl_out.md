# tensorQTL input and output files

This file describes the raw tensorQTL inputs and outputs under:

- `processed-data/11_eQTL_coloc/seurat/`
- `processed-data/11_eQTL_coloc/custom_cluster/`

Each analysis series has its own `tqtl_in/` and `tqtl_out/` directories. Dataset
prefixes are defined by `dataset_id` in each `tqtl_in/prep_manifest.csv`; output
files use that prefix as `<dataset_id>.gene.*`.

## Input files

Each `tqtl_in/` directory contains the expression BED, covariate table, and
expression-PC metadata used by tensorQTL:

- `<dataset_id>.gene.expr.bed.gz`
- `<dataset_id>.gene.covars.txt`
- `<dataset_id>.gene.exprPCs.qs2`
- `prep_manifest.csv`

These files were prepared by `code/11_eQTL_coloc/01_prep_inputs.R`.

Genotype data is read from `processed-data/00_genotypes/plink2/merged_maf05`
by `code/11_eQTL_coloc/02a_tensorQTL_cis.py`.

## tensorQTL output files

Each `tqtl_out/` directory contains raw tensorQTL outputs for its dataset
prefixes:

- `*.gene.cis_qtl_pairs.chr*.parquet`: per-chromosome `map_nominal` outputs.
- `*.gene.map_cis.tab.gz`: permutation-based `map_cis` output, with q-values
  added by the local run script.
- `*.gene.map_independent.txt.gz`: `map_independent` output for features passing
  the script's cis-eQTL FDR threshold.
- `*.gene.DXBPDINT.cis_qtl_top_assoc.txt.gz`: tensorQTL nominal interaction
  top-association output for the available binary BPD diagnosis covariate.
- `*.gene.DXMDDINT.cis_qtl_top_assoc.txt.gz`: tensorQTL nominal interaction
  top-association output for the available binary MDD diagnosis covariate.
- `logs/*.log`: run logs for each dataset prefix.

The interaction files were produced with tensorQTL's nominal interaction mode
for binary Dx-coded covariates available in the covariate table. Treat them as
format-level tensorQTL outputs, not as a final interpretation of a full
three-diagnosis interaction model.

The columns for these tensorQTL formats are documented here:
https://github.com/broadinstitute/tensorqtl/blob/master/docs/outputs.md

## Derived tables

Derived summary tables from these raw tensorQTL outputs are documented in
`processed-data/11_eQTL_coloc/README_tables.md`.
