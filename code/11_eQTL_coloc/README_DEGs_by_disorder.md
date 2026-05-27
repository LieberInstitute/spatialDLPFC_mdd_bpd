# Disorder-Related DEG Lists

This note documents the MDD-related and BPD-related DEG lists used for targeted
eQTL-DEG-GWAS overlap summaries.

## Source Files

The lists are derived from the same standard DEG summary CSVs used by the broad
DEG universe:

- `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

All rows in these final summary files have already passed the omnibus F-test
filter used by this workflow: BH-adjusted `adj.P.Val < 0.05`.

## Definitions

The broad DEG list is the union of all F-test significant genes from the four
source files.

The disorder-related lists further require post-hoc t-test support for the
matching diagnosis-vs-control contrast:

- MDD-related DEG: F-test significant and `F_NTC.MDD_ttest` or
  `M_NTC.MDD_ttest` is significant.
- BPD-related DEG: F-test significant and `F_NTC.BPD_ttest` or
  `M_NTC.BPD_ttest` is significant.

The t-test columns store adjusted p-value bins as text. The workflow treats
`padj<.0001`, `padj<.01`, and `padj<.05` as significant.

## Statistical Rationale

Use F-test only for the broad DEG universe because the omnibus F-test asks
whether a gene has any diagnosis/sex-related signal across the tested contrast
family. This preserves genes with broad, diffuse, or multi-contrast evidence.

Do not use F-test only for MDD-related or BPD-related labels because the F-test
does not identify which disorder contrast drove the signal. A gene can pass the
omnibus F-test because of MDD-vs-NTC, BPD-vs-NTC, MDD-vs-BPD, sex-specific
effects, or combinations.

Use F-test plus contrast t-test for disorder-related lists. This localizes an
already F-test significant gene to MDD-vs-NTC or BPD-vs-NTC evidence.

## Caveats

- F-test only is more sensitive but less specific for disorder labels.
- F-test plus t-test is more specific and better aligned with disorder-specific
  GWAS overlap, but may omit genes where the omnibus F-test is significant and
  no single MDD/BPD post-hoc contrast passes correction.
- These are disorder-related lists, not formal disorder-specific-only lists.
  A gene can be both MDD-related and BPD-related.
