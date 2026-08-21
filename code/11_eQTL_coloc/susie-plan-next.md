# SuSiE next execution plan

## Purpose and target definition

The completed 100-target SuSiE analysis was defined mechanically from the
final `coloc_pass` sheet. It was not a complete list of manuscript-priority or
LD-discussion loci. The 100 comparisons represented 54 genes: 35 BD, 26 MDD,
and 39 SCZD rows. Inclusion required single-causal ABF PP4 > 0.8 at p12 = 1e-5
and H4 > 0.8 for at least half of the plausible p12 grid from 5e-6 through
5e-5. This is a robust-ABF-H4 discovery subset, not a manuscript-priority
target definition.

The next analysis must preserve two separate strata:

1. Confirmatory: the original 100 `coloc_pass` comparisons.
2. Exploratory: a predeclared manuscript/Slack LD-discussion gap set containing
   relevant comparisons for `IFITM2`, `IFITM3`, `SPON2`, and especially
   `SURF1`.

`MAPK3` (five comparisons), `FBLN7` (two), and `PTP4A3` (one) were already in
the completed set. `IFITM2` and `IFITM3` were eQTL examples, `SPON2` was an
exploratory exact-variant GWAS overlap, and `SURF1` was an LD/tie example.
Their different manuscript roles explain their mechanical exclusion but do
not remove the need for locus-specific analysis.

Gene    What the approved no-23andMe ABF analysis showed                 Why excluded
------  ----------------------------------------------------------------  --------------------------------------------------
IFITM2  MDD Astro PP3 = 0.899, PP4 = 0.007; MDD L5 PP3 = 0.986,          Strong evidence for distinct rather than shared
        PP4 = 0.000063                                                    signals
IFITM3  No strong coloc category in any of 24 disorder/context            No H4 result near the inclusion threshold
        comparisons; MDD Astro PP4 = 0.057
SPON2   MDD L6 PP3 = 0.197, PP4 = 0.790; 68% of the plausible prior       Just missed the mandatory primary PP4 > 0.8 rule
        grid passed
SURF1   MDD contexts had PP3 about 0.770 and PP4 below 0.007; no          No shared-signal result; mostly evidence leaning
        context passed                                                    toward distinct signals

The later full-23andMe MDD result for SPON2 L6 (PP4 = 0.928) is out of scope
and must not be used. The approved no-23andMe result was PP4 = 0.790.

## Required pre-execution manifest

Create and review a small manifest before running the exploratory set. For
every row, record the gene, disorder, cell context, GRCh38 locus bounds,
approved GWAS source, and reason for inclusion. At minimum include:

- `IFITM2`: MDD Astro and L5, where coloc.abf favored H3.
- `IFITM3`: manuscript/Slack MDD contexts, including Astro.
- `SPON2`: MDD L6 using only the approved no-23andMe GWAS.
- `SURF1`: all manuscript/Slack-relevant MDD contexts, prioritizing the
  high-H3 comparisons.

Derive the complete context list from the manuscript and Slack record before
viewing new SuSiE results. Do not expand automatically to all 24
disorder/context combinations unless that broader scope is justified.

## Analysis workflow

1. Preserve the reviewed targeted workflow: dense tensorQTL `map_nominal`
   statistics, exact context donors, covariate-adjusted signed eQTL LD, fresh
   approved current no-23andMe GWAS inputs, and ancestry-matched external GWAS
   LD.
2. Apply the same GRCh38 locus definition, allele harmonization, variant
   matching, SuSiE settings, and exclusion rules used for MAPK3 and the
   original targets.
3. Fit eQTL SuSiE directly with genotypes from the relevant donor cohort.
   Reuse an existing eQTL fit only after verifying that its variants, order,
   alleles, statistics, donors, covariates, LD matrix, and settings match.
4. Fit GWAS SuSiE twice against the approved 1000G EUR GRCh38 reference:
   - finite-reference baseline: `R_finite = 503`, `R_mismatch = "none"`;
   - mismatch sensitivity: `R_finite = 503`, `R_mismatch = "eb"`.
5. Run direct susieR and `coloc::runsusie` with identical inputs and settings,
   then run `coloc.susie`. Investigate any fit divergence before interpreting
   colocalization.
6. Preserve all exact genotype ties. Report tied variants as statistically
   indistinguishable rather than as independent evidence.
7. Add locus LD summaries and plots. Report signed r and r-squared from
   unphased donor dosages. Add phased haplotypes and D-prime only when needed
   to answer a population-haplotype question.
8. Keep the exploratory extension and all sensitivity outputs isolated from
   production outputs.

Do not use the full-23andMe GWAS, old coloc caches, or rerun coloc.abf. Existing
coloc.abf results may be used only afterward as descriptive comparators.

## Run order

1. Refit only the MAPK3 GWAS SuSiE pilot with `R_mismatch = "eb"`; retain the
   completed `R_mismatch = "none"` fit as the baseline and reuse verified
   matched eQTL fits.
2. Compare MAPK3 credible sets, PIPs, convergence, mismatch diagnostics,
   runtime, memory, and coloc.susie H3/H4 results between the two GWAS fits.
3. If the MAPK3 fits are reconcilable and diagnostics are acceptable, run both
   GWAS mismatch settings for the predeclared exploratory gap set, with
   `SURF1` prioritized.
4. If those diagnostics remain acceptable, run the EB sensitivity for the
   original 100 comparisons. The existing finite-only results remain the
   comparison baseline.
5. Adapt job parallelism only after measuring per-job memory and CPU use, and
   retain enough memory margin to prevent oversubscription.
6. Perform TOP-LD pairwise validation afterward for decisive SNP pairs and
   credible-set structure. Treat TOP-LD as a separate LD-reference validation,
   not as a replacement for EB mismatch modeling and not as a dense SuSiE LD
   matrix.

## Interpretation and diagnostics

For each exploratory locus, determine whether the original low PP4 or high H3
is best explained by:

- multiple signals violating the locus-wide single-causal ABF assumption;
- genuinely different eQTL and GWAS causal variants;
- insufficient resolution among correlated or tied variants;
- external-LD mismatch detected or accommodated by EB;
- allele, coverage, or reference-panel exclusions;
- weak eQTL or GWAS evidence rather than an LD problem.

Compare donor LD with the ancestry-matched external reference. Inspect
possible recombinant haplotypes and overlay repeats or structural variants if
LD extends unexpectedly far. Treat absence of comparable eQTL and GWAS
credible sets as unresolved, not as negative colocalization evidence.

For each fit, record:

- LD validation and matrix diagnostics;
- allele harmonization, variant counts, and every exclusion with its reason;
- direct susieR and coloc wrapper convergence diagnostics;
- credible sets, lead variants, rsIDs, PIPs, purity, and tied variants;
- `R_mismatch = "none"` versus `"eb"` changes;
- coloc.susie H3/H4 for every eligible signal comparison;
- variants carrying most of the H4 posterior;
- runtime, peak memory, CPU use, package versions, and complete provenance.

Do not describe a result as a shared variant unless H4 supports one causal
variant common to both traits. If posterior support is distributed among
correlated variants, state explicitly that the shared variant's rsID remains
unresolved.

## Deliverables and output separation

Maintain separate directories and summaries for:

1. the original `coloc_pass` results;
2. the EB sensitivity results;
3. the manuscript/Slack exploratory results;
4. the TOP-LD validation.

The final report must compare the SuSiE findings approximately with the
existing approved coloc.abf results for key genes and DEGs. Highlight any
high-H3 ABF cases that SuSiE resolves into separate signals, a supported
shared variant, or an explicitly unresolved result. Report the exploratory
gap set as a manuscript-priority follow-up, not as retrospective members of
the original 100-target analysis.
