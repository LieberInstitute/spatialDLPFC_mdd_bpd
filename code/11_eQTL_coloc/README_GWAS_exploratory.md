# Exploratory GWAS Trifecta Overlaps

This note documents whether SCZD-supported eQTL-DEG-GWASg variants from the
current MBv Seurat eQTL analysis also appear in the exploratory MDD or BPD GWAS
overlap checks.

The companion machine-readable trifecta table is:

- `GWAS_exploratory_trifecta.csv`

The per-variant eQTL overlap tables are written under
`processed-data/11_eQTL_coloc/seurat/tables/`:

- `map_significant_unified_GWASx.csv.gz`
- `map_significant_pairs_GWASx.csv.gz`
- `map_significant_summary_GWASx.csv`

Legacy copies are also written with `GWAS_relaxed` in the filename.

## Thresholds and Context Columns

The SCZD reference uses the strict GWAS threshold, while MDD and BPD matching
uses the most relaxed genotype-matched GWAS cache currently available locally.

| Disorder | Role | Threshold |
|---|---|---:|
| SCZD | reference | `p <= 5e-8` |
| MDD | exploratory matching | `p < 1e-5` |
| BPD | exploratory matching | `p < 1e-5` |

The `p < 1e-5` MDD/BPD cutoff is suggestive/exploratory and relaxed relative
to strict genome-wide significance.

`SCZD_context`, `MDD_context`, and `BPD_context` all refer to MBv Seurat eQTL
contexts from this study, not to a GWAS cohort stratum. Specifically:

- `SCZD_context`: Seurat context(s) where the SCZD-supported eGene/variant is
  part of the broad-DEG overlap.
- `MDD_context` and `BPD_context`: Seurat context(s) where the same variant is
  also observed in the MDD or BPD exploratory broad-DEG exact-variant overlap.
- Blank MDD/BPD columns mean the SCZD variant was not present in the matching
  disorder's `p < 1e-5` genotype-matched GWAS/eQTL-DEG overlap.

In the eQTL tables, strict GWAS columns keep their existing names
(`SCZD_GWAS`, `MDD_GWAS`, `BPD_GWAS`). Exploratory mood-disorder columns use
`GWASx`, for example `MDD_GWASx`, `MDD_GWASxg`, and `MDD_GWASx_p`.

## Simplified Variant Match Table

| Gene | SCZD_context | SCZD_support | SCZD_variant | SCZD_rsid | SCZD_GWAS_p | MDD_context | MDD_variant | MDD_rsid | MDD_GWAS_p | BPD_context | BPD_variant | BPD_rsid | BPD_GWAS_p |
|---|---|---|---|---|---:|---|---|---|---:|---|---|---|---:|
| ARL17B | L2.3 | exact eQTL=SCZD GWAS variant | `chr17:45855941:T:C` | rs7221167 | 1.949979e-08 |  |  |  |  |  |  |  |  |
| ARL17B | Inhb | gene-list-only eQTL variant | `chr17:46025316:C:CT` |  |  |  |  |  |  |  |  |  |  |
| ATF4 | Oligo | gene-list-only eQTL variant | `chr22:39530856:G:A` |  |  |  |  |  |  |  |  |  |  |
| ATF4 | Astro;Inhb;L2.3;L5 | gene-list-only eQTL variant | `chr22:39544222:G:A` |  |  |  |  |  |  |  |  |  |  |
| MAPK3 | Astro;L2.3;L4;L5 | exact eQTL=SCZD GWAS variant | `chr16:30123335:T:C` | rs28529403 | 4.106956e-10 |  |  |  |  | Astro;L2.3;L4;L5 | `chr16:30123335:T:C` | rs28529403 | 7.313917e-08 |
| MAPK3 | Inhb | gene-list-only eQTL variant | `chr16:30311847:G:C` |  |  | Inhb | `chr16:30311847:G:C` | rs148788997 | 5.908946e-07 |  |  |  |  |

## Main Interpretation

The exact SCZD eQTL/GWAS variant overlap that also appears in a mood-disorder
exploratory match is MAPK3 `rs28529403` / `chr16:30123335:T:C`, which appears
in BPD at `p = 7.313917e-08` and in the same MBv contexts as the SCZD row:
Astro, L2.3, L4, and L5.

The MAPK3 Inhb variant `rs148788997` / `chr16:30311847:G:C` appears in MDD at
`p = 5.908946e-07`, but this SCZD row is gene-list-only support rather than an
exact SCZD GWAS variant overlap. ARL17B and ATF4 do not have matching MDD or
BPD variant overlaps in the available `p < 1e-5` caches. Compared with the
prior narrower exploratory run, the wider `p < 1e-5` threshold did not add new
SCZD-supported trifecta variant recoveries.
