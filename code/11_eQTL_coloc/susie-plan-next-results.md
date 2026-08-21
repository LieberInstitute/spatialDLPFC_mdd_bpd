# SuSiE next-plan results

## Scope and authoritative workflow

This run implemented `susie-plan-next.md` without rerunning `coloc.abf`, reading
old coloc caches, using the full-23andMe MDD GWAS, or changing production
outputs. The confirmatory stratum remained the 100 rows mechanically selected
from the final `coloc_pass` table. The exploratory stratum was frozen before
viewing new SuSiE results as nine MDD comparisons: IFITM2 Astro/L5, IFITM3
Astro, SPON2 L6, and SURF1 Astro/Inhb/L4/L5/Oligo.

The workflow was audited against these sources:

- the local `README.md` and `01_prep_inputs.R` for the 119-donor sample order,
  covariates, and dense tensorQTL nominal inputs;
- the official susieR RSS and LD-mismatch documentation:
  https://stephenslab.github.io/susieR/reference/susie_rss.html and
  https://stephenslab.github.io/susieR/articles/rss_mismatch.html;
- the official coloc `runsusie` and `coloc.susie` documentation:
  https://chr1swallace.github.io/coloc/reference/runsusie.html and
  https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html;
- the official 1000 Genomes high-coverage collection and TOP-LD paper/API for
  external LD validation.

Several assumptions required explicit correction. Dense `map_nominal` output
is the correct association-statistic input; the significant and independent
tables are only reporting/target aids. The 119-donor genotypes are the correct
source for covariate-adjusted signed eQTL LD, but they cannot substitute for
ancestry-matched GWAS LD. `R_finite=FALSE` is appropriate for in-sample eQTL
LD, whereas the external 503-EUR GWAS LD uses `R_finite=503` and the EB run adds
`R_mismatch="eb"`. EB can model broad LD mismatch but cannot repair allele
coding errors. Finally, absence of a pure SuSiE credible set is unresolved, not
evidence against colocalization; single-causal ABF does not require a pair of
multi-effect credible sets and can therefore report high PP4 when
`coloc.susie` has no eligible pair.

## Inputs and provenance

- eQTL statistics: full per-gene tensorQTL nominal parquet output.
- eQTL LD: exact context donors, with genotype dosages projected through the
  documented covariate matrix.
- GWAS: current approved BD, SCZD, and MDD no-23andMe BCFs, with SI >= 0.8 and
  exact/swapped allele harmonization.
- GWAS LD: newly extracted target variants from the pinned 1000G high-coverage
  GDS, using 503 unrelated EUR samples. Source GDS MD5:
  `9716fd4fb9562fa2da8edab6a05c4a1e`.
- Model settings: L=10, scaled prior variance 0.2, estimated prior variance,
  Wald z-scores, 95% coverage, minimum absolute credible-set correlation 0.5,
  maximum 1000 iterations, and no refinement.
- Packages: coloc 6.0.1 at commit
  `50fe5291fea7f8ab49823bd86747385d6e56870f`; susieR 0.16.6 at commit
  `ef213feed2cb82419677661a8c986e1504df2c73`.

The exploratory provenance audit hashes 29 direct input files and records all
package, host, source, reference-extraction, harmonization, and exclusion
metadata.

## MAPK3 EB gate

All 24 MAPK3 disorder/context GWAS fits converged under EB. Direct
`susieR::susie_rss` and `coloc::runsusie` were bit-identical for PIPs, alpha,
ELBO, and credible sets. No sensitivity or reliability flag fired.

For the five approved MAPK3 targets, EB preserved every qualitative result:

| Comparison | Baseline H3/H4 | EB H3/H4 | Interpretation |
| --- | --- | --- | --- |
| BD Astro | 0.0787 / 0.9212 | 0.0788 / 0.9212 | shared component; shared SNP unresolved |
| BD L2.3 | 0.1755 / 0.8244 | 0.1756 / 0.8243 | shared component; shared SNP unresolved |
| BD L5 | 0.1267 / 0.8711 | 0.1268 / 0.8710 | shared component; shared SNP unresolved |
| BD Inhb | 0.9855 / 0.0137 | 0.9855 / 0.0137 | distinct eQTL and BD components |
| MDD Inhb | 0.0165 / 0.9825 | 0.0189 / 0.9808 | shared component |

The conditional H4 posterior for the three BD shared rows remains divided
mainly between the dosage-tied rs55732507 and rs28529403 proxies. TOP-LD EUR
independently reports R2=0.998 and D-prime=1.000 for this pair. MDD estimated a
substantial regional mismatch component (effective reference size about 23),
but the MDD Inhb conclusion remained stable and unflagged.

## Exploratory gap set

All 36 direct/wrapper trait fits across nine targets and two mismatch modes
converged and were bit-identical. EB did not change any credible set or
interpretation.

| Gene/context | eQTL CS | GWAS CS | coloc.susie result | Explanation |
| --- | ---: | ---: | --- | --- |
| IFITM2 Astro | 1 | 1 | H3=0.9959, H4=0.000093 | distinct low-LD signals |
| IFITM2 L5 | 1 | 1 | H3=0.9982, H4=0.000037 | distinct low-LD signals |
| IFITM3 Astro | 0 | 1 | no eligible pair | eQTL fit has no retained pure CS |
| SPON2 L6 | 1 | 0 | no eligible pair | GWAS evidence is too diffuse for a retained CS |
| SURF1 Astro | 2 | 0 | no eligible pair | strong eQTL components, no GWAS CS |
| SURF1 Inhb | 1 | 0 | no eligible pair | strong eQTL component, no GWAS CS |
| SURF1 L4 | 2 | 0 | no eligible pair | strong eQTL components, no GWAS CS |
| SURF1 L5 | 2 | 0 | no eligible pair | strong eQTL components, no GWAS CS |
| SURF1 Oligo | 1 | 0 | no eligible pair | broad eQTL CS, no GWAS CS |

IFITM2 directly explains the earlier ABF high-H3 findings: the eQTL lead and
MDD lead are about 64 kb apart but effectively uncorrelated. Reference r2 is
0.000076 for the Astro lead pair and 0.000096 for the L5 lead pair. Thus the
high H3 is signal separation, not an external-LD artifact.

The manuscript IFITM3 variant rs61876236 was explicitly excluded because the
approved MDD input had no matching allele record; it was also absent from the
503-EUR extraction. The remaining matched eQTL fit did not form a pure CS.
This is an allele/coverage limitation plus unresolved eQTL fine-mapping, not a
SuSiE-negative result.

SPON2 rs13119951 leads a 65-variant eQTL CS, but the MDD GWAS fit has no CS.
Its eQTL and GWAS top variants have donor-adjusted r2=0.495 and reference
r2=0.409. TOP-LD EUR reports R2=0.414 and D-prime=0.808. The approved
no-23andMe ABF PP4=0.790 therefore remains suggestive but is not SuSiE
confirmed. The prohibited full-23andMe result was not used.

SURF1 resolves the LD question but not colocalization. The three manuscript
representatives have reference r2=0.919-0.985 and TOP-LD R2=0.910-0.989, so
they form a highly correlated regional structure rather than independent
replications. In contrast, each is nearly uncorrelated with the MDD top signal
rs55924785 about 0.67-0.71 Mb away: reference r2=0.0011-0.0028. EB estimated a
large mismatch component for the SURF1 GWAS region (effective reference size
about 33), but still produced no GWAS CS and no formal reliability flag. This
supports separate regional peaks while leaving multi-signal colocalization
unresolved.

## Original 100-target EB sensitivity

The complete result contains 100/100 targets and 200/200 converged trait fits.
All direct/wrapper comparisons and all recomputed-versus-baseline eQTL reuse
checks are exactly equal. No susieR sensitivity or reliability flag fired.
Baseline had 31 signal pairs in 24 targets; EB had 29 pairs in 23 targets.
GWAS credible-set count changed in five targets; eQTL credible-set count never
changed.

Target-level classifications changed from 76 no-pair, 16 shared, 6 mixed, and
2 distinct at baseline to 77 no-pair, 16 shared, 5 mixed, and 2 distinct under
EB. Most posterior changes were negligible. The important exception was SCZD
L6 HLA-DMA: baseline had one shared pair (H4=0.934) and a second apparent
distinct pair (H3=0.993). EB suppressed the mismatch-sensitive second GWAS
component and left one shared pair (H3=0.0155, H4=0.9845). Its mismatch
estimate reduced the effective reference size to 44 and its maximum penalty
was 5.06, so this is a useful explanation of the high-H3 case but should be
reported as EB-sensitive rather than definitive resolution.

BD L2.3 HCG17 changed from a shared baseline pair (H4=0.935) to no eligible
pair because EB removed its only GWAS CS; its effective reference size was 21.
This result should be downgraded to unresolved. DCC and CPT1C lost extra GWAS
credible sets but already lacked eQTL credible sets, so their no-pair status
did not change.

Across the original table, the strongest remaining high-H3 cases were NAGA
SCZD L6 (H3=0.9886, H4=0.0017) and MAPK3 BD Inhb (H3=0.9855, H4=0.0137), both
supporting distinct components. FADS1 MDD L6, TMEM258 MDD Astro, RNASEH2C SCZD
Astro, and MMD BD Astro/L2.3 retained mixed architectures with both high-H3
and high-H4 signal pairs. These are examples where a locus-wide ABF H4 result
is refined into multiple shared and distinct component comparisons.

Among manuscript DEGs, MAPK3 conclusions were stable. FBLN7 Astro/L2.3 and
PTP4A3 L4 still lacked comparable eQTL/GWAS credible-set pairs and remain
unresolved, not negative or SuSiE-confirmed.

## Diagnostics, exclusions, and resources

Every checkpoint records the signed LD diagnostics, effective rank, finite
reference diagnostics, mismatch lambda, corrected effective reference size,
Q_art, penalties, convergence, warnings, timings, peak RSS, harmonization, and
variant-level exclusions. Exploratory included-variant counts ranged from
2,578 to 4,467. The principal exclusions were no approved GWAS allele match
(548-1,064 per target), absence from the 1000G extraction (8-12), and SI below
0.8 (2-5). No exclusion was hidden.

- MAPK3 EB: 1m55s wall time; 560 MiB peak RSS.
- Exploratory baseline: 1m59s; exploratory EB: 2m55s; peak child RSS 2.23 and
  2.37 GiB, respectively.
- Original 95-row EB batch: 19m24s wall time with adaptive parallelism; 6.40
  total user-CPU hours, median per-process CPU 170%, maximum per-target RSS
  6.41 GiB, and observed aggregate worker RSS about 26 GiB.

TOP-LD returned estimates for 6 of 13 predeclared pairs. The seven absent
pairs are recorded as `not_returned_by_TOP-LD`, not assigned R2=0. The pinned
client SHA256, request parameters, returned rows, and pair audit are retained.
No additional phased donor-haplotype analysis was needed because signed donor
LD, 1000G LD, and TOP-LD R2/D-prime answered the stated pairwise questions.

## Output separation

- MAPK3 EB: `processed-data/11_eQTL_coloc/seurat/coloc/susie_next_no23andMe_20260820/mapk3_eb`
- exploratory baseline/EB and LD:
  `processed-data/11_eQTL_coloc/seurat/coloc/susie_next_no23andMe_20260820/exploratory`
- TOP-LD:
  `processed-data/11_eQTL_coloc/seurat/coloc/susie_next_no23andMe_20260820/topld_validation`
- original 100 EB checkpoints and summaries: `results_eb`, `logs_eb*`, and
  `aggregate_eb_20260820` under the existing `all_targets_20260819` SuSiE run.

Production coloc outputs were not modified. Large checkpoints, compressed
variant ledgers, and plots remain outside git; only code, reviewed manifests,
small summaries, and documentation belong on the `susie` branch.
