# Exact-Tie eQTL Variant Recovery: GWAS, Colocalization, and Trifecta Impact

## Technical summary

The tensorQTL lead/independent tables retained one index variant per reported
signal even when multiple variants had identical donor-level dosages. The
exact-tie recovery audit preserves the original 3,435 signal rows and adds
12,242 previously unreported variant-gene memberships, producing 15,677 rows
across the same signals.

No eGenes were lost or gained. The current and recovered tables both contain
1,390 unique eGenes globally, with zero change in every domain-cell-type
context. The loss affected variant-level connections assigned to existing
eGenes, including exact GWAS overlaps and recognition of strong-coloc lead SNPs
as members of reported eQTL signals.

At the gene-context level, recovery adds 13 exact eQTL-variant/GWAS-variant
overlaps: 5 exploratory MDD, 5 exploratory BD, and 3 strict SCZ. Twelve become
new combined eQTL-GWAS memberships; the thirteenth, SCZ-*ARL17B* in Inhb,
strengthens exact-variant evidence for a membership already present through the
SCZ GWAS gene list. Three new formal eQTL-broad-DEG-GWAS trifecta memberships
are added, all for *SNORC* under exploratory MDD in L4, L5, and Oligo.

Recovery also changes the exact eQTL-signal-member annotation for 11 strong
coloc rows: 10 BD and 1 SCZ, with no MDD changes. These are not new
colocalizations. They are strong colocalizations whose lead SNP was already in
the coloc result but was absent from the one-index eQTL representation.

## Definitions and comparison basis

- **Recovered exact tie:** a variant with the same raw nominal eQTL p-value as
  the reported index and an identical donor-level dosage vector up to allele
  orientation (`abs(r) = 1` within a `1e-12` numerical tolerance).
- **eGene-DEG difecta:** an eGene that is also in the broad DEG union. This is a
  gene-level overlap and is unchanged by variant recovery.
- **eGene-GWAS difecta:** an eGene with combined GWAS support at the stated
  threshold, where combined support is an exact eQTL/GWAS variant match or
  membership in the disorder GWAS gene list. “Difecta” is descriptive wording;
  it is not a stored column name.
- **Formal trifecta:** an eGene in the broad DEG union with combined GWAS
  support. This matches the repository's `trifecta_<DIS>_<level>` definition.
- **Exact GWAS overlap:** an eQTL variant ID that matches a GWAS variant at the
  stated threshold. This excludes gene-list-only support.
- **Strong coloc lead match:** a gated strong-coloc row whose `lead_snp` exactly
  matches a member of the corresponding domain-CT eQTL signal.

The baseline is `map_significant_pairs.csv.gz`. The comparison table is
`map_significant_pairs_tie_recovered.csv.gz`. Both include all lead-cis and
conditionally independent signal rows; no separate cis-only recovery is needed.

## No eGene discovery counts changed

| Domain-CT | Current eGenes | Recovered eGenes | Change |
|---|---:|---:|---:|
| Astro | 561 | 561 | 0 |
| Inhb | 265 | 265 | 0 |
| L2.3 | 925 | 925 | 0 |
| L4 | 211 | 211 | 0 |
| L5 | 770 | 770 | 0 |
| L6 | 461 | 461 | 0 |
| Micro.Vasc | 24 | 24 | 0 |
| Oligo | 79 | 79 | 0 |
| Global unique eGenes | 1,390 | 1,390 | 0 |

The broad-DEG eGene set is likewise unchanged because both the eGene identity
and DEG annotation are gene-level. Therefore no eGene-DEG difecta genes were
lost or rescued. The affected results below arise only when an exact variant ID
is used to connect an existing eQTL signal to another evidence domain.

## Exact GWAS overlaps rescued by disorder and domain-CT

| Disorder and threshold | Domain-CT | Newly recognized exact-overlap genes |
|---|---|---|
| MDD exploratory (`p < 1e-5`) | Inhb | *WNT2B*, *WDR6* |
| MDD exploratory (`p < 1e-5`) | L4 | *SNORC* |
| MDD exploratory (`p < 1e-5`) | L5 | *SNORC* |
| MDD exploratory (`p < 1e-5`) | Oligo | *SNORC* |
| BD exploratory (`p < 1e-5`) | L2.3 | *HLA-DMA*, *HLA-DMB* |
| BD exploratory (`p < 1e-5`) | L5 | *HLA-DMA*, *HLA-DMB*, *AC012213.4* |
| SCZ strict (`p <= 5e-8`) | Astro | *POLR3H* |
| SCZ strict (`p <= 5e-8`) | Inhb | *ARL17B* |
| SCZ strict (`p <= 5e-8`) | L2.3 | *RPS17* |

There are no rescued strict MDD or strict BD exact-overlap memberships.

The union of genes affected in at least one context is:

- MDD exploratory: *WNT2B*, *WDR6*, *SNORC*.
- BD exploratory: *HLA-DMA*, *HLA-DMB*, *AC012213.4*.
- SCZ strict: *POLR3H*, *ARL17B*, *RPS17*.

Genes newly entering the global exact-overlap and combined eGene-GWAS difecta
lists, rather than gaining only an additional context, are:

- MDD exploratory: *WNT2B*.
- BD exploratory: *HLA-DMA*, *HLA-DMB*, *AC012213.4*.
- SCZ strict: *POLR3H*.

*WDR6*, *SNORC*, and *RPS17* already had exact support in another context.
*ARL17B* already had combined SCZ support in Inhb through the GWAS gene list;
recovery adds exact-variant evidence but does not add a new combined-membership
row.

## Formal trifecta changes are limited to SNORC

The only new formal eQTL-broad-DEG-GWAS memberships are:

| Disorder and threshold | Domain-CT | Newly added trifecta gene |
|---|---|---|
| MDD exploratory (`p < 1e-5`) | L4 | *SNORC* |
| MDD exploratory (`p < 1e-5`) | L5 | *SNORC* |
| MDD exploratory (`p < 1e-5`) | Oligo | *SNORC* |

*SNORC* was already a trifecta gene in another context, so no gene is newly
added to the globally collapsed trifecta union. No BD or SCZ formal trifecta
membership changes occur. SCZ-*ARL17B* in Inhb is a broad-DEG exact-overlap
strengthening, but it was already a formal trifecta membership because of its
SCZ GWAS gene-list support.

## Strong-coloc lead-SNP membership rescued by disorder and domain-CT

Exact matching of strong-coloc lead SNPs to reported eQTL signal members rises
from 24 to 35 of the 100 gated strong-coloc rows.

| Disorder | Domain-CT | Strong-coloc genes whose lead SNP is newly recognized as an eQTL signal member |
|---|---|---|
| BD | Astro | *AP001505.1*, *MAPK3* |
| BD | L2.3 | *AP001505.1*, *MAPK3* |
| BD | L5 | *AC012213.4*, *AP001505.1*, *MAPK3*, *RMI2* |
| BD | L6 | *AP001505.1*, *WAC-AS1* |
| SCZ | Inhb | *ARL17B* |

By disorder, BD changes from 8 to 18 exact lead/eQTL-signal matches, MDD remains
3, and SCZ changes from 13 to 14.

The globally affected gene union is:

- BD: *AC012213.4*, *AP001505.1*, *MAPK3*, *RMI2*, *WAC-AS1*.
- SCZ: *ARL17B*.
- MDD: none.

Among these rows, the genes in the broad DEG union are *MAPK3* for BD in Astro,
L2.3, and L5, and *ARL17B* for SCZ in Inhb. This does not create new coloc rows,
and it does not create new formal `trifecta_*` gene memberships; it repairs the
exact lead-SNP-to-eQTL-signal-member annotation.

## Method and robustness checks

The recovery is implemented in [`03ab_recover_indepeQTLs.R`](03ab_recover_indepeQTLs.R).
It runs after `03_eqtl_explore.Rmd`, reads the complete tensorQTL nominal
parquets, and checks dosage equivalence in the correct per-dataset donor set.
It preserves all 3,435 current index rows across all 52 original columns.

Of 15,686 same-p candidate memberships, 15,677 passed the exact-dosage gate.
Nine candidates for one Astro chromosome 2 signal had `r2 = 0.9808` rather than
1 and were excluded. This demonstrates that same p-values or high LD alone were
not accepted as exact ties.

The full audit was run on rsrv16 using staged copies of the current local eQTL
and gated coloc tables because rsrv16's pre-existing processed tables were an
older snapshot. No primary current table was overwritten.

## Interpretation and next use

The one-index tensorQTL tables remain valid summaries of independent eQTL
signals, but they are not exhaustive variant-member lists. They should not be
used directly as exhaustive exact-ID sets for GWAS, coloc-lead, or other
variant-domain intersections.

Use the recovered combined table before eQTL boxplot selection, exact GWAS
overlap summaries, and final coloc lead-membership annotation. The coloc ABF
calculation itself already uses the full nominal parquet data and does not need
to be rerun because of this recovery.

## Audit outputs

- `map_significant_pairs_tie_recovered.csv.gz`: recovered combined eQTL table.
- `map_independent_significant_tie_recovered.csv.gz`: recovered independent-only table.
- `independent_eqtl_tie_recovery_members.csv.gz`: compact signal-member audit.
- `independent_eqtl_tie_recovery_summary_by_context.csv`: recovery scale by context.
- `independent_eqtl_tie_recovery_new_gwas_matches.csv`: rescued exact GWAS signal matches.
- `independent_eqtl_tie_recovery_new_coloc_matches.csv`: rescued strong-coloc lead matches.
- `independent_eqtl_tie_recovery_rejected_same_p_candidates.csv`: excluded non-exact candidates.
