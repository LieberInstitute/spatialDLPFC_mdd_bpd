# GWAS Gene and Variant List Provenance

This file documents the paper-derived gene and variant lists staged beside the BD, MDD, and SCZD GWAS summary statistics. The goal is to make the eQTL/GWAS overlap terminology explicit and to prevent broad method-specific gene lists from being mistaken for high-confidence GWAS gene lists.

Exact variant overlap remains the most precise GWAS evidence. The `GWASg` gene-list overlap is a complementary gene-level bridge for eGenes and DEG overlaps when the eQTL gene is named by a paper-derived high-confidence or prioritized GWAS gene list.

## Recommended Gene Lists for eGene Overlap

| Disorder | Recommended high-confidence/prioritized gene list | Rows | Unique gene symbols | Why this list |
|---|---|---:|---:|---|
| BD | `BD/bpd2024_prioritized_credible_genes.tsv` or identical-source `BD/bpd2024_gene_lists.tsv` | 116 | 116 | Paper-derived credible/prioritized genes from Table S31; compact enough for direct eGene overlap. |
| MDD | `MDD/mdd2025_high_confidence_genes.tsv` | 296 | 295 | Cell Table S8B high-confidence genes; this is the accepted default and the only MDD list that should be treated as high-confidence here. |
| SCZD | `SCZD/sczd2022_prioritized_genes.tsv` | 120 | 120 | Focused Table S12 `Prioritised` list; restrictive, but closest analogue to BD credible/prioritized and MDD high-confidence genes. |

Use these recommended lists when the analysis question is: "Does this eGene match a GWAS high-confidence/prioritized disease gene?" Use exact variant overlap separately when the analysis question is: "Does this eQTL variant itself match a significant GWAS variant?"

## Harmonized Vocabulary

- `variant-only GWAS overlap`: an eQTL variant matches the harmonized significant GWAS variant set for the disorder.
- `prioritized/high-confidence GWAS gene overlap`: an eGene symbol appears in the selected paper-derived high-confidence or prioritized gene list.
- `GWASg`: mixed overlap used in eQTL summaries: variant-only GWAS overlap OR selected GWAS gene-list overlap.
- `broad/source gene lists`: larger paper-derived or method-specific lists useful for provenance and sensitivity analysis, but not necessarily appropriate as the default eGene overlap list.

## How We Decide

For `GWASg`, prefer a paper-curated high-confidence/prioritized gene list over a broad method scan. The gene-list branch of `GWASg` is meant to complement exact variant matching without making the overlap too diffuse.

The selection rule is:

1. Prefer a named high-confidence, prioritized, or credible-gene list from the GWAS paper or its primary supplement.
2. Prefer compact convergent lists over single-method scans with thousands of genes.
3. Keep broad lists documented for provenance and sensitivity checks, but do not call them "high-confidence" unless the paper does.
4. If a disorder has both broad and focused lists, use the focused list as the high-confidence analogue unless the analysis is explicitly a sensitivity analysis.

For SCZD, this means the 120-gene `Prioritised` list is recommended for high-confidence eGene overlap. It is restrictive, but it is the best terminology/method analogue to BD credible/prioritized genes and MDD high-confidence genes. The broader SCZD all-criteria and combined transcriptomic lists remain useful for sensitivity/provenance, not as the default high-confidence list.

## MDD Failure Mode to Avoid

Do not use `MDD/mdd2025_gene_lists.tsv` as the MDD high-confidence list. That file combines fastBAT, Hi-C, and DrugTargetor/MAGMA outputs and has 9,896 rows. It is too broad for a direct high-confidence eGene overlap and was the source of the earlier oversized MDD `GWASg` behavior.

The accepted MDD high-confidence overlap list is:

- source: Cell supplementary `1-s2.0-S0092867424014156-mmc9.xlsx`
- sheet: `Table S8B High-confidence Genes`
- local file: `MDD/mdd2025_high_confidence_genes.tsv`
- current extracted size: 296 rows, 295 unique gene symbols after omitting blank gene values

## Current Code Behavior

This documentation does not change code behavior.

Current constants in `code/11_eQTL_coloc/utils.R` define both source and standardized gene-list roles:

- BD source `broad`: `BD/bpd2024_gene_lists.tsv`
- BD source `prio`: `BD/bpd2024_prioritized_credible_genes.tsv`
- MDD source `broad`: `MDD/mdd2025_high_confidence_genes.tsv`
- MDD source `prio`: `MDD/mdd2025_high_confidence_genes.tsv`
- SCZD source `broad`: `SCZD/sczd2022_gene_lists.tsv`
- SCZD source `prio`: `SCZD/sczd2022_prioritized_genes.tsv`

`loadGWASGeneList(..., use_prio = FALSE)` currently loads the configured `broad` role. For MDD, `broad` and `prio` both point to the high-confidence list. For BD, the `broad` and `prio` source files are effectively the same 116-row credible/prioritized list. For SCZD, `broad` and `prio` differ substantially; this README recommends the 120-gene `prio` file for the high-confidence analogue, but no code default is changed here.

## Shared Extraction Rules

- Gene files use the common schema in `*_gene_lists.tsv` and related focused files.
- Variant files use the common schema in `*_variant_lists.tsv` and related split files.
- Source coordinates were preserved as published. For BD and MDD, source variant coordinates are hg19 where present. For SCZD, prioritized FINEMAP variant coordinates are GRCh37/hg19.
- `hg38_*` fields were filled by matching source rsIDs to the local GRCh38 BCFs with `bcftools query`; allele matching was used when source alleles were available.
- `primary_list=yes` marks a source list that was central/broad in extraction. It does not automatically mean the list is the recommended high-confidence eGene-overlap list.

## Available Gene Lists

### BD 2024/2025

Source paper: O'Connell et al., Nature, doi:10.1038/s41586-024-08468-9.

Source workbook: `41586_2024_8468_MOESM4_ESM.xlsx`
URL: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-08468-9/MediaObjects/41586_2024_8468_MOESM4_ESM.xlsx`

| File | Source table/list | Method label | Rows | Unique gene symbols | Role | Recommendation |
|---|---|---|---:|---:|---|---|
| `BD/bpd2024_gene_lists.tsv` | Table S31, multi-ancestry prioritized credible genes | credible gene prioritization | 116 | 116 | configured source `broad`; same source as focused file | Recommended BD high-confidence/prioritized overlap list. |
| `BD/bpd2024_prioritized_credible_genes.tsv` | Table S31, same rows as above | credible gene prioritization | 116 | 116 | configured source `prio` | Recommended BD high-confidence/prioritized overlap list. |

BD caveats:

- Current exact-variant matching uses the reconstructed integrated European
  `BD/bip2024_eur.hg38.bcf`, including the supplied 23andMe component.
- The integrated BD BCF remains pre-DENTIST and is not the exact final paper
  file. The public `BD/bip2024_eur_no23andMe.hg38.bcf` remains available for
  sensitivity comparisons.
- The paper supplement contains including-self-report/23andMe and
  excluding-self-report/23andMe lists.
- For direct eGene overlap, the 116-row credible/prioritized gene list is appropriately compact; no larger BD gene-list alternative is currently staged.

### MDD 2025

Source paper: Adams et al., Cell, doi:10.1016/j.cell.2024.12.002.

Source GWAS summary Figshare: `https://doi.org/10.6084/m9.figshare.27061255`

Source downstream results Figshare: `https://doi.org/10.6084/m9.figshare.27089614`

| File | Source table/list | Method label | Rows | Unique gene symbols | Role | Recommendation |
|---|---|---|---:|---:|---|---|
| `MDD/mdd2025_high_confidence_genes.tsv` | Cell supplement Table S8B, `High-confidence Genes` | high-confidence gene association | 296 | 295 | configured source `broad` and `prio` | Recommended MDD high-confidence overlap list. |
| `MDD/mdd2025_gene_lists.tsv` | combined fastBAT, Hi-C, and DrugTargetor/MAGMA lists | broad combined source list | 9,896 | 7,681 | provenance/sensitivity | Do not use as high-confidence/default eGene overlap list. |
| `MDD/mdd2025_fastBAT_bonferroni_significant_genes.tsv` | `Online Results (fastBAT).xlsx`, `fastBAT Results` | gene-based association | 1,568 | 1,568 | provenance/sensitivity | Use only for method-specific sensitivity. |
| `MDD/mdd2025_hic_significant_all_tissues_genes.tsv` | `Online Results (hiC).xlsx`, `HiC Gene Associations` | Hi-C gene mapping | 1,034 | 958 | provenance/sensitivity | Use only for method-specific sensitivity. |
| `MDD/mdd2025_drugtargetor_magma_qBH_le_0_05_genes.tsv` | `Online Results (DrugTargetor).xlsx`, `G GENE_results` | MAGMA gene-based association | 7,294 | 6,950 | provenance/sensitivity | Do not use as high-confidence/default eGene overlap list. |

MDD caveats:

- Current exact-variant matching uses the reconstructed integrated European
  `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf`, including the supplied 23andMe
  component. The public no-23andMe BCF remains available for sensitivity
  comparisons.
- The gene and variant lists here were extracted from paper-level downstream Figshare/supplement results, not derived from the local public no23andMe BCF.
- Preserve the source labels (`multi-ancestry`, `European ancestry`, `full_div`, `full_eur`) rather than relabeling downstream lists as with23andMe/no23andMe.

### SCZD 2022

Source paper: Trubetskoy et al., Nature, doi:10.1038/s41586-022-04434-5.

Source supplement zip: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-04434-5/MediaObjects/41586_2022_4434_MOESM11_ESM.zip`

| File | Source table/list | Method label | Rows | Unique gene symbols | Role | Recommendation |
|---|---|---|---:|---:|---|---|
| `SCZD/sczd2022_prioritized_genes.tsv` | Supplementary Table 12, sheet `Prioritised` | focused gene prioritization | 120 | 120 | configured source `prio` | Recommended SCZD high-confidence/prioritized analogue for eGene overlap. |
| `SCZD/sczd2022_gene_lists.tsv` | Table S12 all-criteria genes, Table S12 focused prioritized genes, and SMR/Hi-C support tables | broad combined source list | 1,072 | 703 | configured source `broad` | Use for provenance/sensitivity, not default high-confidence overlap. |
| subset of `SCZD/sczd2022_gene_lists.tsv` | Table S12 all prioritization criteria | gene prioritization criteria | 685 | 682 | source component | Broader alternative if sensitivity is desired. |

SCZD caveats:

- The 120-gene focused prioritized list is intentionally restrictive.
- The broader combined file includes multiple source components: 685 all-prioritization-criteria rows, 120 focused prioritized rows, and SMR/transcriptomic/Hi-C support tables.
- If the analysis goal is maximum sensitivity rather than high-confidence comparability, report the broader list explicitly as `SCZD broad/sensitivity`, not as the default prioritized/high-confidence list.

## Variant List Provenance

Variant files complement gene-list overlaps but serve a different purpose. They are not substitutes for the high-confidence gene-list recommendation above.

| File | Rows | Description |
|---|---:|---|
| `BD/bpd2024_variant_lists.tsv` | 934 | Published locus/signal rows from Tables S5, S6, S8, S11, and S12. |
| `BD/bpd2024_finemapped_credible_variants.tsv` | 295 | Fine-mapped credible variant rows from Tables S28 and S29. |
| `MDD/mdd2025_variant_lists.tsv` | 31,588 | COJO independent signals and fine-mapped credible causal variants. |
| `MDD/mdd2025_cojo_independent_signals.tsv` | 1,319 | COJO selected independent signals from multi-ancestry and European ancestry sheets. |
| `MDD/mdd2025_finemap_credible_causal_variants.tsv` | 30,269 | Fine-mapped credible-causal variants from diverse and EUR analyses. |
| `SCZD/sczd2022_prioritized_variants.tsv` | 1,678 | FINEMAP prioritized nonsynonymous/UTR variants and single-gene credible-set variants from Tables S13 and S14. |

GRCh38 matching status for split variant files:

| File | Match summary |
|---|---|
| `BD/bpd2024_finemapped_credible_variants.tsv` | 289 rsIDs found, alleles not checked; 6 rsIDs not found in local BCF |
| `MDD/mdd2025_cojo_independent_signals.tsv` | 1,311 rsID+allele matches; 8 rsIDs not found in local BCF |
| `MDD/mdd2025_finemap_credible_causal_variants.tsv` | 30,128 rsID+allele matches; 141 rsIDs not found in local BCF |
| `SCZD/sczd2022_prioritized_variants.tsv` | 1,623 rsID+allele matches; 54 rsIDs not found in local BCF; 1 row without rsID |

## Verification Commands

Run these from the staged GWAS directory, for example on `gw`:

```bash
cd ~/work/R/spatialDLPFC_mdd_bpd/processed-data/ref/GWAS
```

Data row counts, excluding the header row:

```bash
for f in \
  BD/bpd2024_gene_lists.tsv \
  BD/bpd2024_prioritized_credible_genes.tsv \
  MDD/mdd2025_high_confidence_genes.tsv \
  MDD/mdd2025_gene_lists.tsv \
  SCZD/sczd2022_gene_lists.tsv \
  SCZD/sczd2022_prioritized_genes.tsv
do
  awk 'END { print FILENAME "\t" NR - 1 }' "$f"
done
```

Unique nonblank `gene_symbol` counts:

```bash
for f in \
  BD/bpd2024_gene_lists.tsv \
  BD/bpd2024_prioritized_credible_genes.tsv \
  MDD/mdd2025_high_confidence_genes.tsv \
  MDD/mdd2025_gene_lists.tsv \
  SCZD/sczd2022_gene_lists.tsv \
  SCZD/sczd2022_prioritized_genes.tsv
do
  printf "%s\t" "$f"
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "gene_symbol") gene_col = i
      next
    }
    gene_col && $gene_col != "" { genes[$gene_col] = 1 }
    END { print length(genes) }
  ' "$f"
done
```

List names, method labels, and `primary_list` values:

```bash
for f in BD/bpd2024_gene_lists.tsv MDD/mdd2025_high_confidence_genes.tsv \
  MDD/mdd2025_gene_lists.tsv SCZD/sczd2022_gene_lists.tsv SCZD/sczd2022_prioritized_genes.tsv
do
  echo "## $f"
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      next
    }
    {
      list[$(h["list_name"])]++
      ev[$(h["evidence_type"])]++
      primary[$(h["primary_list"])]++
      table[$(h["source_table"])]++
    }
    END {
      for (k in list) print "list_name", k, list[k]
      for (k in ev) print "evidence_type", k, ev[k]
      for (k in primary) print "primary_list", k, primary[k]
      for (k in table) print "source_table", k, table[k]
    }
  ' "$f" | sort
done
```

MDD high-confidence sanity check:

```bash
awk 'END {print "MDD high-confidence rows:", NR - 1}' MDD/mdd2025_high_confidence_genes.tsv
awk 'END {print "MDD combined rows:", NR - 1}' MDD/mdd2025_gene_lists.tsv
```

The first command should report 296 rows and the second should report 9,896 rows. If the MDD `GWASg` gene-list branch suddenly behaves like thousands of genes, the combined list was probably used by mistake.
