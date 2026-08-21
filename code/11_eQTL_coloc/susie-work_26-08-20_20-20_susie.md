# SuSiE succession-actions-only work log

## Basis and target definition

The workflow was governed by `susie-plan-next.md`, the local `README.md` and
input-preparation code, current official susieR RSS/mismatch documentation,
current official coloc `runsusie`/`coloc.susie` documentation, the official
1000 Genomes high-coverage resource, and the TOP-LD paper/API. The manuscript
and Slack synthesis defined the exploratory gap before new results were viewed.

The confirmatory targets were the existing 100 `coloc_pass` comparisons: ABF
PP4 > 0.8 at p12=1e-5 and H4 > 0.8 for at least half of the reviewed p12 grid.
The separate exploratory manifest contained nine predeclared MDD comparisons:
IFITM2 Astro/L5, IFITM3 Astro, SPON2 L6, and SURF1
Astro/Inhb/L4/L5/Oligo. It did not retrospectively expand to all contexts.

## Successful actions in execution order

1. Read `susie-plan-next.md`, `README.md`, the manuscript/Slack synthesis, and
   the SURF1 manuscript-selection code. This fixed the target strata and
   verified the all-donor covariate and dense nominal-statistic requirements.

2. Checked current official susieR and coloc documentation. This confirmed
   `R_finite=503` plus `R_mismatch="eb"` for the external GWAS LD sensitivity,
   `R_finite=FALSE` for in-sample eQTL LD, and the direct
   `susie_rss`/`runsusie`/`coloc.susie` sequence.

3. Wrote `susie-exploratory-manifest.tsv` with nine target IDs, gene/context
   pairs, GRCh38 dense-nominal bounds, approved MDD no-23andMe source, and
   evidence-based inclusion reasons. This froze scope before fitting.

4. Added `15_run_MAPK3_gwas_eb.R` and ran it against the 24 verified MAPK3
   pilot objects. The action reused the matched eQTL fits, refit only GWAS with
   EB, compared direct and wrapper fits, ran `coloc.susie`, and wrote isolated
   diagnostics and provenance.

5. Reviewed the MAPK3 gate. All 48 direct/wrapper fits converged and were
   bit-identical; no sensitivity or reliability flag fired, so execution
   advanced under the approved gate.

6. Generalized `08_prepare_targeted_susie_all.R` to accept a reviewed manifest
   and isolated output directory while preserving its 100-target default.
   Preparation successfully extracted nine dense nominal datasets and the
   current approved no-23andMe MDD GWAS records.

7. Reused the checksum-verified 1000G source GDS and ran
   `09_extract_1000g_target_reference.R`. This freshly extracted exact target
   genotypes for 503 unrelated EUR samples on chromosomes 4, 9, and 11 and
   recorded all absent reference variants.

8. Extended `10_run_targeted_susie_one.R` and
   `11_run_targeted_susie_all.sh` with explicit `none`/`eb` modes, isolated
   checkpoint/log directories, direct/wrapper comparisons, EB diagnostics,
   and mismatch metadata.

9. Ran all nine exploratory targets under the finite-reference baseline with
   eight workers. Every job and every trait fit completed and converged; direct
   and wrapper outputs were exactly equal.

10. Ran all nine exploratory targets under EB with eight workers. Every job
    and trait fit completed and converged; no credible-set conclusion changed
    from baseline.

11. Added and ran `16_aggregate_susie_next.R`. It produced reviewable fit,
    credible-set, H4-variant, exclusion, LD, timing, warning, metadata, and
    baseline-versus-EB tables without opening production coloc objects.

12. Added and ran `17_validate_exploratory_ld.R`. It reconstructed raw donor,
    covariate-adjusted donor, and 503-EUR signed LD for reported eQTL members,
    top variants, and credible-set variants, then wrote pair tables, summaries,
    and heatmaps.

13. Added a predeclared TOP-LD pair file and
    `19_run_topld_validation.sh`. The action checksum-pinned the published API
    client, queried EUR LD, retained R2/D-prime results, and explicitly audited
    pairs not returned by the service.

14. Ran `13_record_targeted_susie_provenance.R` for the exploratory directory.
    It hashed 29 direct input files and recorded package commits, host, external
    references, sample panel, harmonization files, and documentation URLs.

15. Started EB sensitivity for the 95 non-MAPK3 rows of the original
    `coloc_pass` set. After measuring about two CPU cores and 1.5-2.2 GiB per
    initial worker, added disjoint small- and middle-locus worker pools to use
    the 32-core host while retaining hundreds of GiB of memory margin.

16. Detected the queue fronts meeting, identified duplicate active target IDs,
    and stopped only the newer copies before checkpoint publication. The older
    processes completed atomically; the main scheduler was resumed and then
    rerun in skip mode to normalize completion state.

17. Verified 95/95 non-MAPK3 EB checkpoints, zero failed final targets, 190/190
    converged non-MAPK3 trait fits, and exact direct/wrapper agreement. Combined
    these with the five completed MAPK3 targets for 100/100 and 200/200 fits.

18. Added and ran `18_aggregate_original_eb.R`. It exported all EB summaries,
    credible sets, conditional H4 variant posteriors, exclusions, diagnostics,
    resources, warnings, and baseline comparisons, and verified exact eQTL
    recomputation versus the reusable baseline fits.

19. Reviewed the ABF comparison read-only. This identified stable distinct
    MAPK3/NAGA signals, stable mixed architectures, the EB-sensitive HLA-DMA
    high-H3 component, and the HCG17 result that became unresolved.

20. Updated `README.md` and wrote `susie-plan-next-results.md`. These document
    commands, output separation, corrected assumptions, provenance, results,
    caveats, resource use, and manuscript-facing interpretation.

21. Parsed every added or modified R script, checked the shell scripts with
    `bash -n`, ran `git diff --check`, and scanned the deliverables for non-ASCII
    characters. All checks passed.

## Final summary

The MAPK3 EB sensitivity converged exactly between direct and wrapper methods
and preserved all five approved conclusions. The exploratory set confirmed
strong H3 for IFITM2; left IFITM3, SPON2, and SURF1 unresolved because one
trait lacked a retained credible set; and explained the SURF1 discussion as a
high-LD eQTL group separated from a distant, low-LD MDD signal. TOP-LD agreed
with the main pairwise LD conclusions.

The original EB sensitivity completed all 100 targets with no nonconvergence,
direct/wrapper divergence, or formal sensitivity/reliability flag. Most
results were stable. HLA-DMA SCZD L6 was the notable high-H3 case explained by
EB: a mismatch-sensitive second GWAS component disappeared, leaving a shared
pair with H4=0.985. HCG17 lost its GWAS credible set and became unresolved.
MAPK3 BD Inhb and NAGA SCZD L6 remained convincing distinct-component cases,
while several loci retained both shared and distinct component pairs.

All new fits, large tables, plots, and logs remain outside git in isolated
SuSiE directories. No production output was changed, no ABF analysis was
rerun, and no prohibited full-23andMe result was used.
