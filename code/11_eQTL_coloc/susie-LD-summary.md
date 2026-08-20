# LD, exact eQTL ties, and targeted multi-signal colocalization

## Scope and sources

This report synthesizes the locus-LD discussion from August 18-20, 2026,
Results 2.6 and Methods 4.8 of the current MBv manuscript, and the completed
targeted SuSiE results. The workflow itself was reviewed against the project
README and tensorQTL preparation, current official coloc and susieR
documentation, and official 1000 Genomes documentation. Full-23andMe GWAS
data and old coloc caches were excluded.

The completed target list was not a manual list of manuscript genes. It was
the 100-row `coloc_pass` sheet from the approved final ABF workbook: 54 genes,
35 BD comparisons, 26 MDD comparisons, and 39 SCZD comparisons. Inclusion
required ABF PP4 > 0.8 at p12 = 1e-5 and H4 > 0.8 across at least half of the
plausible p12 grid from 5e-6 through 5e-5. Consequently, the run covered
robust ABF-H4 findings rather than every manuscript-priority, exact-tie, or
high-H3 locus.

This distinction matters. `MAPK3`, `FBLN7`, and `PTP4A3` were represented in
the 100 targets. `SURF1`, `SPON2`, `IFITM2`, and `IFITM3` were not. The latter
loci require a separate manuscript-priority extension to answer their direct
LD questions; that omission does not mean LD was missing from the 100 runs.

## Was LD calculated for the completed SuSiE analysis?

Yes. Every completed target used a signed, allele-aligned LD matrix matched to
the corresponding summary statistics and variant order.

For eQTL fine mapping, LD came from the exact context-specific cohort: 119
donors in most contexts and 110 in Oligo. It used the same PGEN hard calls as
tensorQTL, mean imputation, and deterministic rank-aware projection through
the documented covariate model. This covariate-adjusted in-sample LD was the
primary matrix; raw donor LD was retained as a sensitivity. Residual variance
was estimated and `R_finite = FALSE` was used because the LD and eQTL
statistics came from the same samples.

For GWAS fine mapping, LD did not come from the 119 donors. It used a fresh
subset of the 2022 1000 Genomes 30x panel: 503 approved unrelated EUR samples,
exact REF/ALT matching, and signed ALT-dosage correlations. Residual variance
was fixed and `R_finite = 503` represented finite external-reference
uncertainty.

Validation included exact allele and variant-order checks, symmetry, unit
diagonal checks, eigenvalue and effective-rank summaries,
`estimate_s_rss()`, `kriging_rss()`, missingness checks, monomorphic-variant
exclusions, and raw-versus-covariate-adjusted donor-LD sensitivity. Twelve
non-MAPK3 eQTL fits with material LD sensitivity received individual-data
diagnostics. Stored and reconstructed eQTL z scores correlated from 0.9873 to
1.0000; rank-aware RSS and correctly compressed individual-level SuSiE
differed by at most 0.0672 in PIP, and interpretable colocalization classes
remained stable.

What was not done was a phased haplotype or D-prime analysis for each locus.
SuSiE-RSS requires signed genotype-predictor correlations, not phased
haplotypes, so this was not a missing SuSiE input. Haplotype mosaics,
recombination-break inspection, and D-prime remain useful follow-ups when the
scientific question is why LD spans an unexpectedly long interval. They are
also still required for the omitted `SURF1`, `SPON2`, `IFITM2`, and `IFITM3`
loci if those specific discussion questions are to be resolved.

## Concrete manuscript loci

Results 2.6 emphasizes a short list rather than the complete eQTL catalog:

- In female MDD astrocyte-dominant spots, `IFITM2` expression is associated
  with `rs10751647` (eQTL q = 0.0072), and `IFITM3` expression is associated
  with `rs61876236` (q = 0.034).
- The colocalization narrative names the mood-disorder DEGs `MAPK3`, `FBLN7`,
  and `PTP4A3`.
- The detailed MAPK3 narrative centers on `rs28529403`, `rs55732507`, and
  `rs148788997`.
- Methods 4.8.3 names an exploratory MDD/SPON2 overlap at `rs13119951`
  (`chr4:1398767:G:T`, MDD p = 9.92e-7).

In the approved GWAS inputs:

- `rs28529403` (`chr16:30123335:T:C`) has SCZ p = 4.11e-10 and BD
  p = 7.31e-8.
- `rs55732507` (`chr16:30130664:T:C`) has SCZ p = 6.43e-10 and BD
  p = 6.97e-8.
- `rs148788997` (`chr16:30311847:G:C`) is the InbCT MAPK3 eQTL highlighted
  for MDD, with eQTL q = 0.00115 and MDD p = 5.91e-7.

The first two SNPs are 7.3 kb apart and have identical observed dosages in the
119-donor cohort. They are strong SCZ variants but only exploratory BD
variants under a strict p < 5e-8 threshold.

## Interpretation of exact ties and long LD

Exact-tie recovery restored variants omitted when tensorQTL reported one
representative from a set with identical 0/1/2 donor dosages. It improves
reporting, but it does not increase independent evidence. For example,
`rs28529403` and `rs55732507` are two labels for one observable MAPK3 genotype
contrast in this cohort. Their recurrence across AstCT, L2/3CT, L4CT, and
L5CT may be biologically informative, but the number of tied SNP labels is not
the number of statistical replications.

Several general points follow:

- Cohort-specific dosage correlation r = 1 does not prove perfect population
  haplotype LD in 1000 Genomes, TOPMed, or the GWAS samples.
- Similar minor-allele frequencies do not establish LD. Equal-frequency minor
  alleles may occur in different donors.
- Matching genotypes at two endpoints do not describe the intervening
  haplotype. Recombinant or switched haplotypes can distinguish candidates.
- Long LD can reflect ordinary low recombination, a small cohort that did not
  sample rare recombinants, selection, repeats, inversions, or other structural
  variation. An informal inspection of the discussed chromosome 9 interval
  found no obvious large indel or inversion, but this does not exclude smaller
  or poorly represented structural features.
- The observed tied SNPs may all tag an untyped, poorly imputed, filtered, or
  structural causal allele. The data establish indistinguishable proxies, not
  that one typed member is necessarily causal.
- SuSiE models multiple effects using association statistics and signed LD. It
  does not require mixed ancestry. Naive ancestry mixing can introduce
  population structure and LD mismatch; ancestry-specific panels and
  stratified sensitivity analyses are preferable.

The `SURF1` example illustrates the unresolved locus question: seven variants
span about 60 kb in AstCT/InbCT, two span about 45 kb in L4CT/L5CT, and six
span about 25 kb in Oligo. Similar MAF across those groups cannot determine
whether they form one LD block. Direct paired dosages are needed for r and
r-squared; phased haplotypes are optional but useful for D-prime,
recombination, and population-LD descriptions.

## MAPK3 multi-signal results

At primary p12 = 1e-5, the five approved MAPK3 targets were:

| Disorder/context | Single-causal ABF result | SuSiE eQTL component | SuSiE GWAS component | Multi-signal interpretation |
|---|---|---|---|---|
| BD AstCT | Lead `rs55732507`; PP4 = 0.920 | 10-member CS, led jointly by tied `rs28529403` and `rs55732507` (PIP = 0.249 each) | 45-member CS led by `rs8054556` (`chr16:29946895:G:A`; PIP = 0.0498), also containing `rs28529403` and `rs55732507` | Shared component pair: H4 = 0.921, H3 = 0.0787 |
| BD L2/3CT | Lead `rs55732507`; PP4 = 0.801 | 16-member CS containing tied `rs28529403` and `rs55732507` (PIP = 0.152 each) | Same 45-member BD CS led by `rs8054556` | Shared component pair: H4 = 0.824, H3 = 0.175 |
| BD L5CT | Lead `rs55732507`; PP4 = 0.827 | 34-member CS containing tied `rs28529403` and `rs55732507` (PIP = 0.173 each) | Same 45-member BD CS led by `rs8054556` | Shared component pair: H4 = 0.871, H3 = 0.127 |
| BD InbCT | Lead `rs55732507`; PP4 = 0.820 | Singleton CS at `rs148788997` (PIP = 0.979) | Same 45-member BD CS led by `rs8054556`, containing `rs28529403` and `rs55732507` | Distinct components: H3 = 0.985, H4 = 0.0137 |
| MDD InbCT | Lead `rs148788997`; PP4 = 0.944 | Singleton CS at `rs148788997` (PIP = 0.978) | 11-member MDD CS led by `rs4787644` (`chr16:30406798:G:A`; PIP = 0.223), with `rs148788997` also present (PIP = 0.0595) | Shared component pair: H4 = 0.982, H3 = 0.0165 |

Here, shared means a coloc-supported regional component pair, not one uniquely
identified shared SNP.

- BD AstCT, L2/3CT, and L5CT pair an eQTL component containing the inseparable
  `rs28529403`/`rs55732507` proxies with a broad BD GWAS component led by
  `rs8054556`. None of those three SNPs is uniquely nominated.
- BD InbCT is different: the eQTL component centers on `rs148788997`, whereas
  the BD GWAS component centers on `rs8054556` and contains the promoter pair.
  SuSiE changes the single-causal H4-favored result to strong H3.
- MDD InbCT supports sharing between the singleton `rs148788997` eQTL
  component and an 11-variant MDD GWAS component led by `rs4787644` that also
  contains `rs148788997`.

The result therefore supports context-dependent MAPK3 signal architecture,
not one shared causal SNP across every context and disorder.

## Other manuscript genes and the target-definition gap

Three other manuscript DEG findings were present in `coloc_pass`:

- MDD/AstCT `FBLN7`, ABF lead `rs13390931`, PP4 = 0.887.
- MDD/L2/3CT `FBLN7`, ABF lead `rs13390931`, PP4 = 0.966.
- SCZ/L4CT `PTP4A3`, ABF lead `rs4129585`, PP4 = 0.947. This was the strict
  same-disorder exact-variant support example.

None produced a comparable pair of SuSiE eQTL and GWAS credible sets. They are
unresolved by this multi-signal analysis, not SuSiE-negative and not
SuSiE-confirmed.

The omitted loci had different ABF roles:

- `IFITM2`: MDD Astro PP3 = 0.899 and PP4 = 0.007; MDD L5 PP3 = 0.986 and
  PP4 = 0.000063. These favor distinct signals rather than sharing.
- `IFITM3`: no strong category across 24 disorder/context comparisons; MDD
  Astro PP4 = 0.057.
- `SPON2`: approved no-23andMe MDD L6 PP3 = 0.197 and PP4 = 0.790. Although
  68% of the plausible prior grid passed, the primary PP4 > 0.8 requirement
  failed. A later full-23andMe result of PP4 = 0.928 is prohibited here and
  must not be substituted.
- `SURF1`: MDD contexts had PP3 near 0.770 and PP4 below 0.007; no context
  passed the H4 gate.

These exclusions explain the 100-row manifest but do not answer the direct LD
questions. The missing work is documented in `TODO.md` as a separate,
manuscript-priority locus extension.

## Possible impacts and explanations

### Statistical impact

- Exact or near-exact collinearity divides PIP among interchangeable variants
  and prevents unique localization even when the regional component is clear.
- A 119-donor cohort may not contain rare recombination events present in the
  population. Wide credible sets can therefore reflect limited LD resolution,
  not many independent causal variants.
- External GWAS LD can mismatch the study population despite correct allele
  harmonization. Finite-reference flags were clear, but mismatch diagnostics
  and a separately reviewed `R_mismatch = "eb"` sensitivity remain advisable.
- Failure to obtain a comparable credible-set pair can arise because one trait
  has no credible set, the regional signal is weak, LD resolution is poor, or
  the traits have different components. It is not automatically evidence
  against colocalization.
- Single-causal ABF can favor H4 when multiple components are present. MAPK3
  BD InbCT demonstrates how multi-signal modeling can instead isolate strong
  H3.

### Biological impact

- Context differences may represent genuinely context-specific regulatory
  components, different expression measurement precision, or different
  effective sample/covariate structures.
- A shared regional component supports compatible genetic architecture but
  does not prove that altered expression mediates disease risk.
- A tied set can identify a regulatory haplotype without identifying the
  causal nucleotide or molecular mechanism.
- Structural variants, repeats, or an untyped allele may explain apparently
  excessive SNP-level multiplicity and should be investigated when ordinary
  dosage LD is insufficient.

### Reporting impact

- Counts of variants, genes, contexts, and disorder comparisons require
  separate denominators. Several tied SNPs in several contexts do not equal
  several independent loci.
- Exact strict SCZ support for `rs28529403` and `rs55732507` must not be
  described as strict BD support. Their highlighted BD colocalizations use
  exploratory BD associations.
- Tied variants should be reported as one statistically inseparable signal
  group with all member IDs preserved.
- ABF and SuSiE results should be described as colocalization evidence, not as
  proof of a causal variant or a causal expression-to-disease relationship.

## Recommended wording

> At the MAPK3 locus, LD-aware multi-signal colocalization supported shared
> BD-GWAS and MAPK3-eQTL components in AstCT, L2/3CT, and L5CT. The eQTL
> credible sets retained the dosage-identical promoter-proximal variants
> rs28529403 and rs55732507 and could not distinguish between them. In InbCT,
> however, the MAPK3 eQTL component centered on rs148788997 was distinct from
> the BD-GWAS component. The rs148788997-centered InbCT eQTL component was
> compatible with the MDD-GWAS signal.

> Exact-tie recovery improves variant reporting but does not create
> independent evidence. Variants with identical cohort dosages should be
> presented as statistically indistinguishable signal members. Multi-signal
> colocalization can separate regional components only where the association
> statistics and LD contain discriminatory information.

## Further work

1. Create a second target manifest for `SURF1`, `SPON2`, `IFITM2`, and
   `IFITM3`, limited to manuscript- and discussion-relevant contexts and
   disorders. Keep it separate from the completed ABF-H4 target set.
2. Apply the same approved dense-summary, allele-harmonization, signed-LD, and
   SuSiE settings. Do not use significance-filtered variants, full-23andMe
   inputs, or old coloc caches.
3. Produce locus LD summaries using donor dosage r and r-squared, including
   exact-tie membership, context sample size, missingness, allele order, and
   comparison with ancestry-matched external LD.
4. Where long LD remains unexplained, inspect phased haplotypes, recombinant
   donors, D-prime, recombination maps, repeats, and structural-variant
   annotations. Keep these population-genetic diagnostics distinct from the
   signed matrix required by SuSiE-RSS.
5. Run prior sensitivity and a separately reviewed external-LD mismatch
   sensitivity. Use ancestry-specific panels rather than a naively mixed
   reference.
6. Update manuscript language and supplemental results to distinguish shared,
   distinct, mixed, and unresolved component results. Avoid unique causal-SNP
   claims for tied or broad credible sets.

## Bottom line

The completed 100-target run did calculate the LD required by SuSiE correctly.
It established that the `rs28529403`/`rs55732507` MAPK3 eQTL component shares
regional support with a BD GWAS component in AstCT, L2/3CT, and L5CT; that the
`rs148788997` InbCT eQTL is distinct from the BD component; and that the same
`rs148788997` eQTL component is compatible with the MDD signal.

The remaining gap is target definition. The robust-ABF-H4 `coloc_pass` set was
not equivalent to the manuscript/LD-discussion locus set. Direct analyses of
`SURF1`, `SPON2`, `IFITM2`, and `IFITM3` remain necessary before their LD
questions can be considered resolved.
