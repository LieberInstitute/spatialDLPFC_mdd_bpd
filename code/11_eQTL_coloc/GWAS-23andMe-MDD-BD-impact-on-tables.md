# Impact of 23andMe-inclusive MDD and BD GWAS on eQTL/DEG tables

Generated: 2026-07-29 01:57:16 EDT

## Scope

This comparison holds the 119-donor genotype data, tensorQTL eQTL results, DEG definitions, curated GWAS gene lists, coloc priors, and sensitivity gate fixed. It changes only the MDD and BD GWAS inputs from public European no-23andMe statistics to reconstructed European statistics that include 23andMe.

The deltas therefore measure the effect of switching supplied GWAS files, not an isolated marginal effect of the 23andMe cohorts. The integrated meta-analysis statistics and available row-level QC fields also differ from the public files; in particular, integrated MDD/BD has no combined imputation-quality (`SI`) field.

## Terminology

- `no23` means the archived public European GWAS that excludes 23andMe. `full` means the local 23andMe-inclusive European reconstruction obtained by meta-analyzing that public result with the delivered 23andMe component; it does not mean an official final paper file. MDD is a release-matched reconstruction. BD uses v7.0 associations with v7.2 annotations and is `pre-DENTIST`, meaning the paper's final post-meta-analysis LD-based DENTIST QC was not reproduced.
- A `canonical variant` is identified as GRCh38 `chromosome:position:REF:ALT`. An `eQTL row` or `eQTL pair` is one tensorQTL association between such a variant and an eGene, meaning a gene whose expression is associated with that variant, in a cell context and analysis split. The eQTL result, not the GWAS, supplies the gene assignment.
- `strict` or `exact_variant` means that the same canonical eQTL variant has GWAS P <= 5e-8. `suggestive_p1e5` means exact-variant GWAS P < 1e-5.
- `broad_DEG` is membership in the project-wide union of the author-recommended DEG tables. `MDD_DEG` and `BD_DEG` are the diagnosis-specific NTC-versus-MDD and NTC-versus-BD DEG subsets.
- A `curated GWAS gene` comes from a paper-derived GWAS gene list held fixed between runs. `variant_or_curated_gene` means exact variant support OR curated-gene-list membership, so it is not necessarily a variant-level overlap.

## GWAS inputs

| disorder | version | release | bcf | sha256 | SI_filter |
| --- | --- | --- | --- | --- | --- |
| MDD | no23andMe | public_no23andMe | pgc-mdd2025_no23andMe_eur_v3-49-24-11.hg38.bcf | 4fd59d1a7cdb5cc47f7c24aa0ca0fb377249432e8616d3a2e8157fd796f6bdcd | >=0.8 |
| MDD | full23andMe | full23andMe | pgc-mdd2025_eur_v3-49-24-11.hg38.bcf | a8b30df37b032097920ded06697aa60851f0078e681f1a595b6d5f3ec70347f3 | none |
| BD | no23andMe | public_no23andMe | bip2024_eur_no23andMe.hg38.bcf | cb4e246c5a570d3dad74327ff5c4df50f0edd7b9915ac0f4ce59a397b403adaf | >=0.8 |
| BD | full23andMe | full23andMe_preDENTIST | bip2024_eur.hg38.bcf | 1d502351659d81aa503da64779b7cc9e80c6819bca9e2ab5d5ac6efc3bb531be | none |

The integrated MDD file is the 23andMe-inclusive European meta-analysis. The integrated BD file is the reconstructed 23andMe-inclusive European meta-analysis and remains pre-DENTIST; it is not represented here as the exact final paper release. SI is absent from both integrated BCFs, so no post-integration SI filter was applied. The archived public files used SI >= 0.8.

## Cache and table validation

| table | old_rows | new_rows | row_delta | stable_eqtl_deg_columns_equal | numeric_tolerance |
| --- | --- | --- | --- | --- | --- |
| map_significant_pairs | 3,435 | 3,435 | 0 | yes | 1e-12 |
| nominal_BH05 | 455,143 | 455,143 | 0 | yes | 1e-12 |
| map_GWASx | 3,435 | 3,435 | 0 | yes | 1e-12 |

Stable non-GWAS eQTL and DEG columns were equal across all compared 03-series row tables at the documented 1e-12 numeric tolerance: yes.

Curated broad and prioritized GWAS gene-list files were byte-identical between runs:

| disorder | list | n_genes_no23 | n_genes_full | sha256_identical |
| --- | --- | --- | --- | --- |
| BD | gene_list | 116 | 116 | TRUE |
| BD | prio_gene_list | 116 | 116 | TRUE |
| MDD | gene_list | 295 | 295 | TRUE |
| MDD | prio_gene_list | 295 | 295 | TRUE |

Matched variants in the fixed 119-donor PLINK2 target:

| disorder | p1e5_no23 | p1e5_full | p1e5_delta | p5e8_no23 | p5e8_full | p5e8_delta |
| --- | --- | --- | --- | --- | --- | --- |
| BD | 12,660 | 35,186 | +22,526 | 2,894 | 10,470 | +7,576 |
| MDD | 44,555 | 93,770 | +49,215 | 14,053 | 35,668 | +21,615 |

These are matched variant rows, not independent GWAS loci.

## eQTL/GWAS and DEG-scoped overlap changes

### Exact eQTL-GWAS variant overlaps (no DEG filter)

This table counts eQTL rows whose exact canonical variant passes the stated GWAS threshold; `rows_*` are context-specific eQTL association rows and `genes_*` are distinct tensorQTL eGenes, with no DEG or curated-GWAS-gene requirement. `map_significant_pairs` is the significant cis/independent eQTL set, whereas `nominal_BH05` is the broader nominal eQTL set passing BH FDR 0.05.

| table | disorder | level | rows_no23 | rows_full | row_delta | genes_no23 | genes_full | gene_delta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| map_significant_pairs | BD | strict | 6 | 11 | +5 | 4 | 4 | 0 |
| map_significant_pairs | BD | suggestive_p1e5 | 12 | 24 | +12 | 7 | 11 | +4 |
| map_significant_pairs | MDD | strict | 8 | 26 | +18 | 6 | 12 | +6 |
| map_significant_pairs | MDD | suggestive_p1e5 | 24 | 67 | +43 | 15 | 39 | +24 |
| nominal_BH05 | BD | strict | 752 | 2,147 | +1,395 | 13 | 29 | +16 |
| nominal_BH05 | MDD | strict | 2,593 | 6,466 | +3,873 | 36 | 86 | +50 |

### DEG-scoped significant-eQTL overlaps

This table restricts significant eQTL rows by DEG status: `exact_variant` with `broad_DEG` is the exact GWAS-variant/eQTL/DEG trifecta, while `variant_or_curated_gene` with a diagnosis-specific DEG permits either exact-variant support or curated-GWAS-gene support. `gene_names_*` lists the distinct eGenes in each intersection.

| disorder | level | annotation | scope | rows_no23 | rows_full | row_delta | gene_names_no23 | gene_names_full |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BD | strict | exact_variant | broad_DEG | 0 | 0 | 0 |  |  |
| BD | strict | variant_or_curated_gene | BD_DEG | 0 | 0 | 0 |  |  |
| BD | suggestive_p1e5 | exact_variant | broad_DEG | 4 | 0 | -4 | MAPK3 |  |
| BD | suggestive_p1e5 | variant_or_curated_gene | BD_DEG | 4 | 0 | -4 | MAPK3 |  |
| MDD | strict | exact_variant | broad_DEG | 0 | 4 | +4 |  | SNORC |
| MDD | strict | variant_or_curated_gene | MDD_DEG | 0 | 0 | 0 |  |  |
| MDD | suggestive_p1e5 | exact_variant | broad_DEG | 6 | 6 | 0 | MAPK3, SNORC, SPON2 | MAPK3, SNORC, SPON2 |
| MDD | suggestive_p1e5 | variant_or_curated_gene | MDD_DEG | 1 | 1 | 0 | SPON2 | SPON2 |

Exact gained/lost row and gene lists are in `processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28/03_exact_variant_annotation_changes.tsv`. Combined variant-or-curated-gene and DEG-scope counts are in `03_overlap_metrics_comparison.tsv` in the same directory.

## Colocalization changes

| disorder | candidates_no23 | candidates_full | candidate_delta | raw_strong_no23 | raw_strong_full | gated_strong_no23 | gated_strong_full | gated_delta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BD | 584 | 679 | +95 | 35 | 22 | 35 | 22 | -13 |
| MDD | 769 | 1,280 | +511 | 28 | 46 | 26 | 46 | +20 |

Cell contexts with a changed sensitivity-gated strong-coloc count:

| disorder | context | gated_no23 | gated_full | delta |
| --- | --- | --- | --- | --- |
| BD | Astro | 7 | 2 | -5 |
| BD | Inhb | 4 | 3 | -1 |
| BD | L4 | 1 | 0 | -1 |
| BD | L5 | 7 | 4 | -3 |
| BD | L6 | 5 | 2 | -3 |
| MDD | Astro | 3 | 4 | +1 |
| MDD | Inhb | 3 | 6 | +3 |
| MDD | L2.3 | 7 | 11 | +4 |
| MDD | L5 | 5 | 12 | +7 |
| MDD | L6 | 2 | 7 | +5 |

Sensitivity-gated strong-coloc gains and losses:

| disorder | status | loci | genes |
| --- | --- | --- | --- |
| BD | gained | 17 | AL049840.5, AL645608.7, ASPDH, CACNA1B, FADS1, HLA-DMA, LINC01954, LRRC37A2, MED24, PPFIA1, TMEM106B, TMEM258 |
| BD | lost | 30 | AC012213.4, AP001505.1, CCS, CDHR1, HCG17, MAPK3, MMD, NTSR1, RMI2, RPRD2, SERPINI1, VWA5B2, WAC-AS1 |
| MDD | gained | 28 | AL049840.5, AL596257.1, CKS2, IQCB1, MAEL, MGLL, MGMT, MYOM2, PPFIA1, PTP4A3, RBM23, RETREG2, SCLY, SLC25A12, SLC2A11, SPON2, TSFM, TUBGCP6 |
| MDD | lost | 8 | AC010857.1, FBLN7, MAP3K7, PNMA8A, SPATA20, SPSB2 |

Posterior agreement for candidate loci present in both runs:

| disorder | shared_candidates | PP4_pearson | median_abs_PP4_change | max_abs_PP4_change | lead_SNP_unchanged_pct |
| --- | --- | --- | --- | --- | --- |
| BD | 312 | 0.766 | 0.000 | 0.886 | 25.3 |
| MDD | 628 | 0.897 | 0.000 | 0.798 | 47.1 |

Candidate-row totals are screening categories, not counts of independent GWAS loci. The sensitivity-gated strong calls are the appropriate primary comparison for the final coloc tables.

## Reproducibility

Frozen archive: `/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/archive/no23andMe_2026-07-28`

Machine-readable comparison tables: `/home/gpertea/work/R/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28`

Regenerate this report after rerunning 03-05 with:

```bash
Rscript code/11_eQTL_coloc/07_compare_23andme_impact.R
```
