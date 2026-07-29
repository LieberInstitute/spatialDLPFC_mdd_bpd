# Exploratory GWAS Overlaps

This note documents how exploratory MDD and BD GWAS exact-variant overlaps are
represented in the current Seurat eQTL tables.

Current MDD and BD exact-variant statistics use the integrated European BCFs
containing the supplied 23andMe component. The BD reconstruction remains
pre-DENTIST. Paper-derived curated gene lists are unchanged by this GWAS
summary-statistics update.

The authoritative downstream tables are under
`processed-data/11_eQTL_coloc/seurat/tables/`:

- `map_significant_summary.csv`
- `map_significant_pairs.csv.gz`
- `map_significant_unified.csv.gz`

The companion SCZD-focused example table in this directory is:

- `GWAS_exploratory_trifecta.csv`

## Thresholds

| Disorder | Strict overlap | Exploratory overlap |
|---|---:|---:|
| SCZD | `p <= 5e-8` | not used |
| MDD | `p <= 5e-8` | `p < 1e-5` |
| BD | `p <= 5e-8` | `p < 1e-5` |

The `p < 1e-5` MDD/BD cutoff is an exploratory exact-variant check. It is
integrated into the same row-level and summary tables as the strict overlaps.

## Table Columns

Row-level tables use the cleaned GWAS overlap convention:

- `DIS_gwasVar_strict`: exact eQTL variant overlap with the strict GWAS set.
- `DIS_gwasVar_exp`: exact eQTL variant overlap with the exploratory GWAS set;
  present for MDD and BD.
- `DIS_gwasGene`: curated GWAS gene-list overlap.
- `DIS_gwas_strict`: `DIS_gwasVar_strict OR DIS_gwasGene`.
- `DIS_gwas_exp`: `DIS_gwasVar_exp OR DIS_gwasGene`; present for MDD and BD.
- `DIS_gwasP`, `DIS_gwasBeta`, `DIS_gwasBetaSE`: exploratory exact-variant
  GWAS statistics for MDD/BD variant matches.

`map_significant_summary.csv` reports combined eGene counts:

- `n_DIS_gwas_strict`: unique eGenes with strict variant or gene-list support.
- `n_DIS_gwas_exp`: unique eGenes with exploratory variant or gene-list support;
  present for MDD and BD.

Variant-only counts are not separate summary columns; use the row-level
`DIS_gwasVar_*` flags when exact variant support is needed.

## Context Columns In The Example Table

`SCZD_context`, `MDD_context`, and `BD_context` all refer to Seurat eQTL
contexts from this study, not to GWAS cohort strata.

- `SCZD_context`: Seurat context(s) where the SCZD-supported eGene/variant is
  part of the broad-DEG overlap.
- `MDD_context` and `BD_context`: Seurat context(s) where the same variant is
  also observed in the MDD or BD exploratory broad-DEG exact-variant overlap.
- Blank MDD/BD columns mean the SCZD variant was not present in the matching
  disorder's `p < 1e-5` genotype-matched GWAS/eQTL-DEG overlap.

## Simplified Variant Match Table

| Gene | SCZD_context | SCZD_support | SCZD_variant | SCZD_rsid | SCZD_gwasP | MDD_context | MDD_variant | MDD_rsid | MDD_gwasP | BD_context | BD_variant | BD_rsid | BD_gwasP |
|---|---|---|---|---|---:|---|---|---|---:|---|---|---|---:|
| ARL17B | L2.3 | exact eQTL=SCZD GWAS variant | `chr17:45855941:T:C` | rs7221167 | 1.949979e-08 |  |  |  |  |  |  |  |  |
| ARL17B | Inhb | gene-list-only eQTL variant | `chr17:46025316:C:CT` |  |  |  |  |  |  |  |  |  |  |
| ATF4 | Oligo | gene-list-only eQTL variant | `chr22:39530856:G:A` |  |  |  |  |  |  |  |  |  |  |
| ATF4 | Astro;Inhb;L2.3;L5 | gene-list-only eQTL variant | `chr22:39544222:G:A` |  |  |  |  |  |  |  |  |  |  |
| MAPK3 | Astro;L2.3;L4;L5 | exact eQTL=SCZD GWAS variant | `chr16:30123335:T:C` | rs28529403 | 4.106956e-10 |  |  |  |  |  |  |  |  |
| MAPK3 | Inhb | gene-list-only eQTL variant | `chr16:30311847:G:C` |  |  | Inhb | `chr16:30311847:G:C` | rs148788997 | 1.825744e-06 |  |  |  |  |

## Main Interpretation

No exact SCZD eQTL/GWAS variant in this simplified table also has an integrated
MDD or BD same-variant match at `p < 1e-5`. The prior no-23andMe BD match for
MAPK3 `rs28529403` / `chr16:30123335:T:C` is absent from the integrated BD
cache.

The MAPK3 Inhb variant `rs148788997` / `chr16:30311847:G:C` appears in the
integrated MDD results at `p = 1.825744e-06`, but that SCZD row is gene-list
support rather than an exact SCZD variant overlap. ARL17B and ATF4 do not have
matching integrated MDD or BD variant overlaps in the `p < 1e-5` caches.
