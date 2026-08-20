# Targeted SuSiE workflow review

Reviewed 2026-08-19 against the installed source and current official pages:

- https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
- https://chr1swallace.github.io/coloc/reference/runsusie.html
- https://chr1swallace.github.io/coloc/reference/coloc.susie.html
- https://stephenslab.github.io/susieR/reference/susie_rss.html
- https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html
- https://stephenslab.github.io/susieR/articles/finemapping.html
- https://internationalgenome.org/data-portal/data-collections/1000genomes_30x/
- https://www.internationalgenome.org/faq/which-datasets-include-related-individuals/

## Corrections and decisions

1. `coloc::runsusie()` documentation still says `prior_variance`, but installed
   `susieR::susie_rss()` 0.16.6 accepts `scaled_prior_variance`. Passing
   `prior_variance` is an error, not a supported partial match. The pilot uses
   `scaled_prior_variance = 0.2` as the initial value with
   `estimate_prior_variance = TRUE`.
2. The primary fits use z scores, so quantitative-trait `sdY` scaling and the
   case-control `s` parameter do not enter SuSiE. `runsusie()` requires one
   scalar `N`; eQTL uses the actual donor count and each GWAS uses the median
   variant-specific `NE`, with its range and variation reported as an
   approximation forced by the coloc interface.
3. Both eQTL and GWAS z scores are Wald statistics, so `z_method = "wald"`.
   `L = 10`, 95% coverage, minimum absolute CS correlation 0.5, estimated
   prior variance, and no refinement are explicit. Residual variance is
   estimated for covariate-adjusted in-sample eQTL LD (`R_finite = FALSE`), as
   recommended by current susieR documentation. It remains fixed for external
   GWAS LD. A missing 95% credible set is retained as a result.
4. Official susieR guidance recommends covariate-adjusted in-sample LD when
   covariates were removed from the univariate regressions. tensorQTL removes
   covariates from both genotype and phenotype, so raw donor LD is not the
   matched primary matrix.
5. The prepared covariate files contain an explicit intercept. tensorQTL
   1.0.10 centers every covariate column and then uses an unpivoted PyTorch QR.
   The centered intercept is zero, but PyTorch still returns a full Q column.
   The extra null-space vector is implementation- and device-dependent and is
   not part of the intended regression model. Primary eQTL LD therefore uses a
   deterministic, rank-aware projection on the complete covariate design. This
   is equivalent to retaining an intercept in the model while omitting its
   zero centered column. Exact tensorQTL-style CPU QR and raw donor LD are
   declared sensitivities. Stored tensorQTL summary statistics remain primary
   and are checked against rank-aware and tensorQTL-style reconstructed
   statistics; existing tensorQTL outputs are not regenerated.
6. The tensorQTL wrapper loaded hard calls with `PgenReader.load_genotypes()`,
   not imputed dosages. Primary eQTL LD therefore uses the same hard calls and
   mean imputation. Using `pgenlibr::Read()` dosages would be mismatched.
7. GWAS LD cannot use the 119 brain donors. The candidate reference contains
   the 503 EUR samples from the original 2,504 unrelated phase-3 panel, drawn
   from the 3,202-sample GRCh38 30x call set. Genotype ALT counts provide
   signed dosage correlations; haplotype phase is not used.
8. The current susieR development documentation supports finite-reference
   uncertainty through `R_finite`. Primary fits use `R_finite = 503` for the
   external GWAS panel and `R_finite = FALSE` for exact in-sample eQTL LD.
   Because `coloc.susie(..., susie.args=)` cannot pass different arguments to
   its two traits, the dataset-input route is only an implementation check
   under common trusted-R/fixed-residual settings. Primary colocalization uses
   separately fitted trait-appropriate objects. `R_mismatch` remains `"none"`;
   it will not be used to conceal failed harmonization or poor LD diagnostics.
   The 2026-07-28 official mismatch vignette now recommends considering both
   `R_finite` and `R_mismatch = "eb"` for routine external-reference analyses.
   This pilot and its same-method target extension retain the reviewed
   finite-only primary model; the EB model is a separately reviewed future
   sensitivity, not an automatic correction. All primary finite-reference
   sensitivity and reliability flags must be reported.
9. `SI >= 0.8` is enforced strictly. Missing SI is excluded, despite the
   older helper accepting missing SI. There is no GWAS p-value filter.
10. Exact and near-exact genotype ties are not pruned or collapsed. They are
    logged with a reversible membership table. No `nearPD`, shrinkage,
    arbitrary pruning, or reduced credible-set coverage is used.
11. `susieR::estimate_s_rss()` and `kriging_rss()` are diagnostics without an
    official universal failure threshold. They are reported with outliers and
    plots; allele/order checks and direct tensorQTL correlation checks remain
    mandatory.
12. `runsusie()` is a wrapper around `susie_rss()` plus annotation, not an
    independent fine-mapping algorithm. The direct-versus-wrapper comparison
    tests input and argument translation. If equivalent, later targets need
    only `runsusie()` because `coloc.susie()` needs its annotation.

No GPU backend is documented for `coloc.susie()` or `susie_rss()`. Scaling is
by independent CPU processes across reusable region/trait fits, with BLAS and
data.table threads capped per process.
