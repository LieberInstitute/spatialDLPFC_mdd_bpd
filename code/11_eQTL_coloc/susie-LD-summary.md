# LD, exact eQTL ties, and targeted multi-signal colocalization

## Scope and sources

This report synthesizes the locus-LD discussion from August 18-20, 2026,
Results 2.6 and Methods 4.8 of the current MBv manuscript, the completed
targeted SuSiE analysis, the empirical-Bayes (EB) LD-mismatch sensitivity, and
the manuscript-priority exploratory extension. The workflow was reviewed
against the project README and tensorQTL preparation, current official coloc
and susieR documentation, and official 1000 Genomes and TOP-LD documentation.
Full-23andMe GWAS data and old coloc caches were excluded, and existing
`coloc.abf` results were used only as read-only quasi-validation.
Detailed execution results and package/input provenance are recorded in
`susie-plan-next-results.md` and `susie-work_26-08-20_20-20_susie.md`.

The completed target list was not a manual list of manuscript genes. It was
the 100-row `coloc_pass` sheet from the approved final ABF workbook: 54 genes,
35 BD comparisons, 26 MDD comparisons, and 39 SCZD comparisons. Inclusion
required ABF PP4 > 0.8 at p12 = 1e-5 and H4 > 0.8 across at least half of the
plausible p12 grid from 5e-6 through 5e-5. Consequently, the run covered
robust ABF-H4 findings rather than every manuscript-priority, exact-tie, or
high-H3 locus.

This distinction matters. `MAPK3`, `FBLN7`, and `PTP4A3` were represented in
the 100 confirmatory targets. `SURF1`, `SPON2`, `IFITM2`, and `IFITM3` were not
in that ABF-H4-selected set. A separate manifest was therefore frozen before
viewing new SuSiE results and analyzed nine MDD comparisons: IFITM2 AstCT and
L5CT, IFITM3 AstCT, SPON2 L6CT, and SURF1 AstCT, InbCT, L4CT, L5CT, and Oligo.
These exploratory results do not retrospectively change the confirmatory
target definition.

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
uncertainty. All 100 confirmatory targets and all nine exploratory targets
were also fitted with `R_mismatch = "eb"` as a separately reported
sensitivity.

Validation included exact allele and variant-order checks, symmetry, unit
diagonal checks, eigenvalue and effective-rank summaries,
`estimate_s_rss()`, `kriging_rss()`, missingness checks, monomorphic-variant
exclusions, and raw-versus-covariate-adjusted donor-LD sensitivity. Twelve
non-MAPK3 eQTL fits with material LD sensitivity received individual-data
diagnostics. Stored and reconstructed eQTL z scores correlated from 0.9873 to
1.0000; rank-aware RSS and correctly compressed individual-level SuSiE
differed by at most 0.0672 in PIP, and interpretable colocalization classes
remained stable.

The exploratory loci additionally received raw donor, covariate-adjusted
donor, and 503-EUR pairwise LD summaries. Thirteen predeclared pairs were
submitted to TOP-LD: six returned estimates and seven were explicitly recorded
as not returned, never as zero LD. TOP-LD agreed with the principal MAPK3,
SPON2, and SURF1 conclusions.

What was not done was a complete phased-haplotype or D-prime analysis for
every locus. SuSiE-RSS requires signed genotype-predictor correlations, not
phased haplotypes, so this was not a missing model input. Haplotype mosaics,
recombination-break inspection, and structural-variant analysis remain useful
follow-ups when the scientific question is why LD spans an unexpectedly long
interval.

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

The completed `SURF1` validation answers the original LD question. The three
manuscript representative variants are one highly correlated regional
structure: 503-EUR r-squared = 0.919-0.985 and TOP-LD EUR R2 = 0.910-0.989.
They are not independent replications. Each is nearly uncorrelated with the
distant MDD lead `rs55924785`: 503-EUR r-squared = 0.0011-0.0028. Phased
haplotypes remain optional for describing recombination and D-prime but are
not needed to establish this separation.

## MAPK3 multi-signal results

At primary p12 = 1e-5, the five approved MAPK3 targets were:

| Disorder/context | Single-causal ABF result | SuSiE eQTL component | SuSiE GWAS component | Multi-signal interpretation |
|---|---|---|---|---|
| BD AstCT | Lead `rs55732507`; PP4 = 0.920 | 10-member CS, led jointly by tied `rs28529403` and `rs55732507` (PIP = 0.249 each) | 45-member CS led by `rs8054556` (`chr16:29946895:G:A`; PIP = 0.0498), also containing `rs28529403` and `rs55732507` | EB H4 supports one shared causal variant, but its identity is unresolved: H4 = 0.9212, H3 = 0.0788 |
| BD L2/3CT | Lead `rs55732507`; PP4 = 0.801 | 16-member CS containing tied `rs28529403` and `rs55732507` (PIP = 0.152 each) | Same 45-member BD CS led by `rs8054556` | EB H4 supports one shared causal variant, but its identity is unresolved: H4 = 0.8243, H3 = 0.1756 |
| BD L5CT | Lead `rs55732507`; PP4 = 0.827 | 34-member CS containing tied `rs28529403` and `rs55732507` (PIP = 0.173 each) | Same 45-member BD CS led by `rs8054556` | EB H4 supports one shared causal variant, but its identity is unresolved: H4 = 0.8710, H3 = 0.1268 |
| BD InbCT | Lead `rs55732507`; PP4 = 0.820 | Singleton CS at `rs148788997` (PIP = 0.979) | Same 45-member BD CS led by `rs8054556`, containing `rs28529403` and `rs55732507` | EB-distinct components: H3 = 0.9855, H4 = 0.0137 |
| MDD InbCT | Lead `rs148788997`; PP4 = 0.944 | Singleton CS at `rs148788997` (PIP = 0.978) | 11-member MDD CS led by `rs4787644` (`chr16:30406798:G:A`; PIP = 0.223), with `rs148788997` also present (PIP = 0.0595) | EB H4 supports one shared causal variant, but its identity is unresolved: H4 = 0.9808, H3 = 0.0189 |

Here, shared has the standard H4 meaning: the eQTL and GWAS signals are
inferred to have one causal variant in common. H4 does not necessarily identify
which SNP that is. Its posterior support can be distributed across several
correlated credible-set variants, including variants that are statistically
indistinguishable in this cohort.

- In BD AstCT, L2/3CT, and L5CT, H4 supports one causal variant shared by the
  MAPK3 eQTL and BD GWAS signals. The eQTL credible sets contain the
  inseparable `rs28529403`/`rs55732507` proxies, while the broad BD GWAS
  credible set is led by `rs8054556`. None of those three SNPs is uniquely
  nominated as the shared variant.
- BD InbCT is different: the eQTL component centers on `rs148788997`, whereas
  the BD GWAS component centers on `rs8054556` and contains both promoter SNPs.
  SuSiE changes the single-causal H4-favored result to strong H3.
- MDD InbCT supports sharing between the singleton `rs148788997` eQTL
  component and an 11-variant MDD GWAS component led by `rs4787644` that also
  contains `rs148788997`.

The result therefore supports context-dependent MAPK3 signal architecture,
not one shared causal SNP across every context and disorder.

## EB sensitivity across the 100 confirmatory targets

All 200 EB trait fits converged, and each direct susieR result was numerically
identical to its coloc wrapper result. Baseline analysis yielded 31
credible-set pairs in 24 targets; EB yielded 29 pairs in 23 targets. eQTL
credible-set counts did not change. GWAS credible-set counts changed in five
targets, and no formal finite-reference or mismatch reliability flag fired.

The baseline classes were 76 no-pair, 16 shared, six mixed, and two distinct;
the EB classes were 77 no-pair, 16 shared, five mixed, and two distinct. Most
conclusions were stable. The main explanatory change was HLA-DMA SCZ L6CT:
EB removed a mismatch-sensitive second GWAS component, eliminating the
baseline high-H3 pair and leaving a shared pair with H4 = 0.9845. HCG17 BD
L2/3CT instead lost its only GWAS credible set and was downgraded from H4 =
0.935 to unresolved. NAGA SCZ L6CT and MAPK3 BD InbCT remained convincing
distinct-component results.

## Other manuscript genes and exploratory extension

Three other manuscript DEG findings were present in `coloc_pass`:

- MDD/AstCT `FBLN7`, ABF lead `rs13390931`, PP4 = 0.887.
- MDD/L2/3CT `FBLN7`, ABF lead `rs13390931`, PP4 = 0.966.
- SCZ/L4CT `PTP4A3`, ABF lead `rs4129585`, PP4 = 0.947. This was the strict
  same-disorder exact-variant support example.

None had both an eQTL credible set and a GWAS credible set that `coloc.susie`
could compare. This is compatible with high single-causal ABF PP4 because
`coloc.abf` integrates regional single-variant evidence without requiring each
trait to produce a retained, pure SuSiE credible set. These findings are
unresolved by multi-signal analysis, not SuSiE-negative and not
SuSiE-confirmed.

The completed exploratory extension gave the following results:

- `IFITM2`: AstCT yielded H3 = 0.9959 and H4 = 0.000093; L5CT yielded H3 =
  0.9982 and H4 = 0.000037. The respective eQTL and GWAS leads are effectively
  uncorrelated, confirming distinct signals.
- `IFITM3`: AstCT had one GWAS credible set but no retained eQTL credible set.
  The manuscript variant `rs61876236` lacked an approved MDD GWAS allele match
  and was absent from the 503-EUR extraction, so the result remains unresolved.
- `SPON2`: L6CT had one eQTL credible set but no retained GWAS credible set.
  Donor-adjusted r-squared = 0.495, 503-EUR r-squared = 0.409, and TOP-LD EUR
  R2 = 0.414 between the eQTL lead `rs13119951` and MDD variant `rs6851528`.
  Approved no-23andMe ABF PP4 = 0.790 remains suggestive but is not SuSiE
  confirmed.
- `SURF1`: AstCT, InbCT, L4CT, L5CT, and Oligo all produced one or two eQTL
  credible sets but no retained 95% GWAS credible set. EB did not change this.

For SURF1, the approved MDD GWAS itself is large and the top regional variant
`rs55924785` has imputation quality 0.978, but its p value is 3.0e-7 and its
GWAS PIP is only 0.177 without EB and 0.184 with EB. The unfiltered 95%
component contains 94 variants whose minimum absolute correlation is 0.00014,
so it correctly fails the default purity filter. An exploratory 80% set
contains six highly correlated variants with minimum absolute correlation =
0.987, showing a coherent core plus a diffuse posterior tail. This is limited
localization, not evidence that the GWAS file is generally low quality.

SURF1 also has meaningful external-LD mismatch: EB estimated a corrected
effective reference size near 33 rather than 503. Nevertheless, EB did not
recover a credible set, and independent LD sources confirm separation from
the eQTL components. The current evidence therefore leans strongly toward
distinct eQTL and MDD signals while remaining formally unresolved by
`coloc.susie`.

## Possible impacts and explanations

### Statistical impact

- Exact or near-exact collinearity divides PIP among interchangeable variants
  and prevents unique localization even when the regional component is clear.
- A 119-donor cohort may not contain rare recombination events present in the
  population. Wide credible sets can therefore reflect limited LD resolution,
  not many independent causal variants.
- External GWAS LD can mismatch the study population despite correct allele
  harmonization. The completed EB sensitivity showed that most results were
  stable, while HLA-DMA and HCG17 required revised interpretation.
- Failure to obtain comparable eQTL and GWAS credible sets can arise because
  one trait has no credible set, the regional signal is weak, LD resolution is
  poor, or the traits have different components. It is not automatically
  evidence against colocalization.
- Single-causal ABF can favor H4 when multiple components are present. MAPK3
  BD InbCT demonstrates how multi-signal modeling can instead isolate strong
  H3.

### Biological impact

- Context differences may represent genuinely context-specific regulatory
  components, different expression measurement precision, or different
  effective sample/covariate structures.
- H4 support is compatible with one variant driving both associations, but it
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
> colocalization can distinguish multiple causal signals only where the
> association statistics and LD contain discriminatory information.

> At SURF1, the manuscript eQTL representatives form one highly correlated
> regional structure but are nearly uncorrelated with the distant leading MDD
> association. The approved MDD signal did not yield a retained 95% pure SuSiE
> credible set, so multi-signal colocalization remains unresolved rather than
> negative. Existing ABF and LD evidence favor distinct regional signals.

## Further work

1. For SURF1, seek LD matched to the contributing MDD GWAS cohorts or a larger
   ancestry-matched reference. Actual GWAS-cohort LD would be preferable.
2. Treat 80% and 90% SURF1 credible-set analyses only as transparent
   sensitivity checks; do not replace the predeclared 95% primary result.
3. Revisit SURF1 when a larger approved non-23andMe MDD GWAS becomes available
   and test whether it localizes the distant signal or reveals a secondary
   signal in the eQTL block.
4. Replicate the SURF1 eQTL components in a larger independent brain or
   cell-type-specific cohort. Phasing and recombinant-donor inspection may
   then help distinguish the tied eQTL proxies.
5. Update manuscript language and supplemental results to distinguish shared,
   distinct, mixed, and unresolved component results. Avoid unique causal-SNP
   claims for tied or broad credible sets.

## Bottom line

The completed confirmatory and exploratory runs used the LD required by
SuSiE, converged identically between direct susieR and coloc wrappers, and
were stable under EB mismatch correction in most loci. MAPK3 retained shared
BD components in AstCT, L2/3CT, and L5CT, a distinct BD InbCT architecture,
and a shared MDD InbCT component. HLA-DMA showed how EB can remove a
mismatch-sensitive apparent H3 component, while HCG17 became unresolved.

The manuscript-priority extension confirmed distinct IFITM2 signals, left
IFITM3 and SPON2 unresolved because one trait lacked a retained credible set,
and resolved the SURF1 LD structure without establishing colocalization.
SURF1's eQTL representatives form a high-LD group separated from a distant,
weakly localized MDD association. Better GWAS-matched LD and stronger
association information are needed for a definitive multi-signal result; the
present evidence favors H3 over H4.
