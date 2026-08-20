# SuSiE locus-extension TODO

The completed 100-target SuSiE analysis was defined by the final `coloc_pass`
sheet. It was not a complete list of manuscript-priority or LD-discussion loci.
The following omitted loci require the same approved LD-aware treatment if the
goal is to answer those locus-specific questions.

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

## Required work

1. Define a manuscript/LD-discussion target manifest separately from the
   ABF-H4-selected `coloc_pass` manifest. Include at least `SURF1`, `SPON2`,
   `IFITM2`, and `IFITM3`; `MAPK3`, `FBLN7`, and `PTP4A3` are already covered.
2. Select the relevant disorder/context comparisons from the manuscript and
   locus discussion before execution. Do not expand automatically to all 24
   disorder/context combinations unless that broader scope is justified.
3. Use the reviewed targeted workflow: dense tensorQTL `map_nominal`
   statistics, exact context donors, covariate-adjusted signed eQTL LD, fresh
   approved no-23andMe/current GWAS inputs, and ancestry-matched external GWAS
   LD. Do not use full-23andMe inputs or old coloc caches.
4. Preserve all exact genotype ties. Report tied variants as statistically
   indistinguishable rather than as independent evidence.
5. Add direct locus LD summaries and plots. Report signed r and r-squared from
   unphased donor dosages; add phased haplotypes and D-prime only if the
   population-haplotype question requires them.
6. Compare donor LD with an ancestry-matched external panel, inspect possible
   recombinant haplotypes, and overlay repeats or structural variants when LD
   extends unexpectedly far.
7. Treat absence of an eQTL/GWAS credible-set pair as unresolved, not as a
   negative colocalization result.
8. Keep this extension isolated from production outputs and report it as a
   manuscript-priority follow-up, not as part of the completed 100-target run.

## Scope correction

The 100 completed comparisons represented 54 genes: 35 BD, 26 MDD, and 39
SCZD rows. Inclusion required single-causal ABF PP4 > 0.8 at p12 = 1e-5 and
H4 > 0.8 for at least half of the plausible p12 grid from 5e-6 through 5e-5.
This is a robust-ABF-H4 discovery subset, not a manuscript-priority target
definition.

The manuscript colocalization genes already represented were `MAPK3` (five
comparisons), `FBLN7` (two), and `PTP4A3` (one). `IFITM2` and `IFITM3` were
eQTL examples, `SPON2` was an exploratory exact-variant GWAS overlap, and
`SURF1` was an LD/tie example. Their different roles explain the mechanical
exclusion but do not remove the need for a locus-specific extension.

The later full-23andMe MDD result for SPON2 L6 (PP4 = 0.928) is out of scope
and must not be used. The approved no-23andMe result was PP4 = 0.790.
