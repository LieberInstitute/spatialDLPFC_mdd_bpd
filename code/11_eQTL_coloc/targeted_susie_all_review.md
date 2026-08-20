# Targeted SuSiE colocalization review

Date: 2026-08-20

## Scope and completion

- Completed all 100 rows in the final `coloc_pass` table: 5 reviewed MAPK3
  pilot rows and 95 new checkpointed runs.
- All 200 primary trait fits converged. No worker failed.
- No `coloc.abf()` model was rerun. Existing ABF columns were read only for
  quasi-validation.
- No full-23andMe GWAS, old coloc cache, or production output was read or
  modified. BD and MDD used the approved no23andMe BCFs; SCZD used the approved
  current public BCF.

## Reviewed method

The eQTL fit used stored dense tensorQTL `map_nominal` beta/SE statistics and
signed LD from the exact context donors' PGEN hard calls. Genotypes were mean
imputed and projected through the rank-aware documented covariate model.
Residual variance was estimated because this is in-sample LD. The 119 donor
genotypes were not used for GWAS LD.

GWAS LD used a fresh subset of the 2022 1000 Genomes 30x phased panel: 503
approved unrelated EUR samples, exact REF/ALT matching, ALT allele counts,
fixed residual variance, and `R_finite = 503`. The source GDS was freshly
downloaded from the University of Washington SeqArray catalog's conversion of
the official panel. Its MD5 is `9716fd4fb9562fa2da8edab6a05c4a1e`.

The current susieR documentation now also recommends considering
`R_mismatch = "eb"` for broader external-reference population mismatch. It was
not made primary here because the reviewed MAPK3 pilot explicitly fixed
`R_mismatch = "none"`, all subsequent targets were authorized to use the same
method, and this option must not substitute for allele QC. Every finite-panel
reliability flag was false. The elevated external-LD mismatch diagnostics are
still reported below; an EB mismatch run would be a separately reviewed
sensitivity analysis, not a silent replacement.

MAPK3 established that direct `susieR::susie_rss()` and `coloc::runsusie()`
were bit-identical under matched inputs and arguments for all 48 pilot trait
fits. Therefore the 95 later runs used the faster `runsusie()` path and passed
trait-specific prefit objects to `coloc.susie()`.

## Primary results

Of 100 target comparisons:

| SuSiE interpretation at p12 = 1e-5 | Targets |
|---|---:|
| no eQTL/GWAS credible-set pair | 76 |
| at least one strong shared pair (H4 >= 0.8), no strong distinct pair | 16 |
| both strong shared and strong distinct pairs | 6 |
| strong distinct pair only (H3 >= 0.8) | 2 |

There were 31 credible-set pair comparisons across 24 targets. Overall, 22
targets had at least one H4 >= 0.8 pair and 8 had at least one H3 >= 0.8 pair.
The 6 mixed cases were:

| Disorder | Context | Gene | max H4 | max H3 |
|---|---|---|---:|---:|
| BD | Astro | MMD | 0.9830 | 0.9859 |
| BD | L2.3 | MMD | 0.9914 | 0.9711 |
| MDD | Astro | TMEM258 | 0.9672 | 0.9840 |
| MDD | L6 | FADS1 | 0.9525 | 0.9921 |
| SCZD | Astro | RNASEH2C | 0.9892 | 0.9838 |
| SCZD | L6 | HLA-DMA | 0.9338 | 0.9927 |

These are useful examples where single-causal ABF H4 summarized the region as
shared, while SuSiE found one shared pairing and another pairing consistent
with distinct signals.

## Read-only ABF quasi-validation and high-H3 cases

This is not independent validation: ABF and SuSiE reuse the same association
statistics. It tests whether the single-causal ABF pattern is compatible with
multi-signal fine mapping.

Among the 20 final-table rows with ABF PP3 >= 0.10, SuSiE classified 9 as a
strong shared pair, 2 as a strong distinct pair, and 9 as having no pair of
credible sets to compare.

Two cases are especially informative:

- BD Inhb MAPK3: ABF PP4 = 0.820 and PP3 = 0.162, but SuSiE gives H3 =
  0.9855 and H4 = 0.0137. This changes the interpretation to distinct eQTL and
  GWAS signals after multiple-signal fine mapping.
- SCZD L6 NAGA: ABF PP4 = 0.878 and PP3 = 0.103, but SuSiE reports two strong
  distinct-signal pairings (max H3 = 0.9886) and max H4 = 0.00139.

Other higher-H3 ABF rows were often clarified in the opposite direction. BD
L2.3 MAPK3 retained shared support (SuSiE H4 = 0.824), as did CORO7, RPRD2,
HCG17, RPS17, and BD Astro TMEM258. Rows with no pair should not be called
resolved: one or both traits lacked a usable credible set under the reviewed
model.

For DEG targets, 4 of 8 rows had strong shared support, 1 had strong distinct
support, and 3 had no pair. MAPK3 was shared in BD Astro, BD L2.3, BD L5, and
MDD Inhb, but distinct in BD Inhb. FBLN7 in MDD Astro/L2.3 and PTP4A3 in SCZD
L4 had no pair.

## Prior and LD sensitivity

The cross-trait prior matters. Strong H4 targets numbered 22 at p12 = 1e-5,
22 at p12 = 1e-4, but only 5 at p12 = 1e-6. Strong H3 targets remained 8 at
all three p12 values. H4 magnitudes should therefore be reported with the p12
sensitivity table, not as prior-free probabilities.

Rank-adjusted versus raw-donor eQTL LD produced identical credible sets for 84
of 100 rows. All 12 non-MAPK3 cases with a credible-set change or PIP shift
greater than 0.1 received an individual-data residual-space diagnostic. The
rank/raw LD matrices differed by as much as 0.364, explaining large raw-LD
PIP shifts. Stored versus rank-reconstructed z scores remained highly
correlated (0.9873 to 1.0000). Rank-aware RSS and correctly compressed
individual SuSiE differed by at most 0.0672 in PIP.

The primary colocalization interpretation was stable across primary, raw-LD,
rank-reconstructed, and compressed-individual routes for the sensitivity
cases that yielded interpretable pairs (BD Astro TMEM258, MDD Astro TMEM258,
MDD Oligo OTUD7A, and SCZD Astro RNASEH2C). Other sensitive cases had no
eQTL/GWAS credible-set pair, so their PIP change did not generate a positive
colocalization claim.

## LD validation, exclusions, and warnings

- Fresh 1000 Genomes extraction: 223,294 of 224,411 requested allele-specific
  variants matched (99.50%); 1,117 were absent. All 503 samples matched the
  approved order and no reference genotypes were missing.
- Comparison-level ledger totals: 364,746 included; 57,576 no GWAS allele
  match; 4,199 SI below 0.8; 741 absent from the exact 1000 Genomes subset.
  There were no additional monomorphic exclusions.
- eQTL `estimate_s_rss`: median 0.000145, range 0.0000235 to 0.00131.
- GWAS `estimate_s_rss`: median 0.0321, range 0.00688 to 0.1528, consistent
  with the expected greater uncertainty of external LD. The 6 regions above
  0.08 all had no coloc credible-set pair.
- All GWAS `R_finite` sensitivity and reliability flags were false. ASPDH and
  MYOM2 were the only fits just above the documentation's 0.2
  effective-rank/reference-size guideline (0.204 and 0.208); both had no pair.
- Recorded model warnings were either the expected coloc no-credible-set
  message or coloc's low-minimum-p input check. Beta variance was supplied as
  SE squared, so the latter is informational for weak regions.

## Runtime, memory, and provenance

The 95-run batch used 8 checkpointed workers on the 32-CPU host. Wall time was
about 9 minutes; summed worker user CPU was 1.62 hours. Median peak RSS was
1.61 GiB and maximum peak RSS was 5.75 GiB. The MAPK3 pilot separately took
20 minutes and peaked at 3.70 GiB. The fresh GDS reference extraction took
3 minutes 56 seconds and peaked at 2.65 GiB.

Primary package builds were coloc 6.0.1 and susieR 0.16.6 from the pinned
MAPK3 library. The audit directory contains hashes for all 137 direct input
files, package paths/commits, host/session information, source URLs, and the
reference extraction manifest.

## Conclusions

The complete targeted run supports 22 regions with at least one strong shared
signal pair, but it also materially refines the ABF interpretation: 6 regions
contain both shared and distinct pairings, and MAPK3 BD Inhb plus NAGA SCZD L6
favor distinct signals. Most ABF-pass rows (76/100) did not yield a pair of
credible sets and should not be carried forward as SuSiE-confirmed coloc.
Positive H4 calls are sensitive to the p12 prior and should be presented with
that sensitivity and the external-LD diagnostics.
