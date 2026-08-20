# SuSiE methods review

Reviewed 2026-08-19 for the targeted MAPK3 `coloc.susie` pilot.

## Which SuSiE method is most widely used?

The two R entry points under consideration are not competing fine-mapping
methods:

- `susieR::susie_rss()` is the canonical SuSiE implementation for regression
  summary statistics.
- `coloc::runsusie()` is a wrapper around `susieR::susie_rss()`. It repeats
  fitting until convergence and adds variant labels and metadata used by
  `coloc.susie()`.
- `coloc::coloc.susie()` is a downstream signal-specific colocalization method,
  not a separate SuSiE fitter. It compares log Bayes factors from the SuSiE
  fits for two traits.

Consequently, `susieR` is the dominant reference implementation for SuSiE
fine-mapping and diagnostics. In coloc workflows, `runsusie()` is the usual
convenience interface because its result is ready for `coloc.susie()`.

Running direct `susie_rss()` and `runsusie()` on identical effective inputs and
settings is an implementation and argument-translation check. Their PIPs,
credible sets, ELBO, log Bayes factors, convergence state, and variant order
should agree numerically. If that equivalence is established in the MAPK3
pilot, later coloc targets can use `runsusie()` alone while retaining the
direct fit as a pilot validation record.

Sources:

- coloc `runsusie()` reference:
  <https://chr1swallace.github.io/coloc/reference/runsusie.html>
- coloc `coloc.susie()` reference:
  <https://chr1swallace.github.io/coloc/reference/coloc.susie.html>
- coloc SuSiE vignette:
  <https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html>
- canonical susieR project and citations:
  <https://stephenslab.github.io/susieR/>

## Python and GPU implementations

Broad's current tensorQTL repository contains relevant GPU-enabled SuSiE
functionality:

- The command-line program exposes `cis_susie` and `trans_susie` modes.
- The SuSiE implementation is written in Python using PyTorch.
- It selects CUDA automatically when a compatible GPU is available.
- It operates on individual-level genotype, phenotype, and covariate data.

This is scientifically relevant to an eQTL workflow because the exact donor
genotypes, expression phenotypes, and covariates are available. It is not,
however, a drop-in replacement for the current `susieR::susie_rss()` workflow:

- It is an individual-level implementation, not the complete current
  summary-statistics RSS interface.
- It does not expose the complete set of current susieR finite-reference and
  LD-mismatch diagnostics.
- Its output is not directly accepted by `coloc.susie()`.
- Exact numerical and feature parity with current susieR should not be assumed
  without a dedicated comparison.

Recent pure-Python SuSiE ports also exist, including NumPy/SciPy and optional
Numba implementations. They are new and are not maintained by the original
Stephens Lab susieR team. They should not replace susieR as the primary
publication implementation without extensive independent validation.

For a MAPK3-sized region, GPU acceleration is unlikely to be decisive. Allele
harmonization, genotype extraction, LD construction, and diagnostics may take
more time than the SuSiE fit itself. GPU tensorQTL becomes more attractive when
many phenotypes or loci are fine-mapped jointly from individual-level data.
CPU parallelism across independent regions remains the simpler established
approach for the targeted coloc workflow.

Sources:

- Broad tensorQTL repository:
  <https://github.com/broadinstitute/tensorqtl>
- tensorQTL command modes, including `cis_susie` and `trans_susie`:
  <https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/tensorqtl.py>
- tensorQTL PyTorch SuSiE implementation:
  <https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/susie.py>
- Recent pure-Python port example:
  <https://github.com/omicverse/py-susie>

## Recommended SuSiE practices

### 1. Prefer exact individual-level data when available

For the eQTL trait, `susieR::susie()` on the exact donor genotype matrix and
expression phenotype is the strongest reference fit. It avoids external LD
reference mismatch and provides an independent check of the RSS reconstruction.

For the targeted pilot, compare:

1. Direct individual-level `susie()`.
2. Direct `susie_rss()` using matched covariate-adjusted in-sample LD.
3. `runsusie()` using exactly the same RSS inputs and settings.

Agreement among these fits is more informative than substituting an
unvalidated Python or GPU port.

Source:

- susieR individual-level fine-mapping vignette:
  <https://stephenslab.github.io/susieR/articles/finemapping.html>

### 2. Use signed and allele-aligned LD

SuSiE requires the signed correlation matrix `r`, not `r2`, absolute
correlation, or an unsigned LD measure. Variant names and order must be
identical between the summary statistics and LD matrix.

Require explicit agreement on:

- genome build and chromosome coordinates;
- variant identity and ordering;
- REF and ALT alleles;
- effect allele and beta or z-score direction;
- multiallelic representation;
- missing, monomorphic, and excluded variants.

Do not infer allele orientation from position alone. Palindromic alleles need
frequency-supported resolution or exclusion when orientation is ambiguous.

Source:

- susieR summary-statistics fine-mapping vignette:
  <https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html>

### 3. Prefer matched in-sample LD

The preferred LD matrix is computed from the exact genotype matrix and samples
used to obtain the association statistics. If the association regression
removed covariate effects from genotypes, susieR recommends computing the
in-sample LD from genotypes after those same effects have been removed.

For the eQTL analysis this means:

- Use the exact context-specific donors.
- Use the same genotype representation used by tensorQTL association testing.
- Apply the same missing-genotype treatment.
- Residualize genotypes against the intended context-specific covariate model.
- Calculate signed correlations from the residualized genotype matrix.

The covariate design must be full rank. A redundant explicit intercept should
be represented statistically but should not be allowed to create an arbitrary
null-space QR column.

The 119 brain donors are appropriate for eQTL LD but are not an appropriate
LD reference for an external GWAS. GWAS LD should be study matched when
legitimately available, or otherwise ancestry matched and explicitly treated
as reference-panel LD.

Source:

- susieR summary-statistics guidance, including covariate-adjusted in-sample
  LD:
  <https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html>

### 4. Retain dense regional summary statistics

Fine-mapping is a joint regional model. Include all variants in the defined
candidate region that pass predeclared genotype and imputation QC. Do not
select variants using GWAS or eQTL significance before fitting SuSiE.

Define the region before examining SuSiE output and ensure that each candidate
signal has sufficient surrounding sequence and LD context. Record the region
definition and every exclusion.

### 5. Treat sample size and LD construction as model inputs

Record and validate:

- association sample size;
- context-specific eQTL donor count;
- GWAS effective sample-size convention;
- LD reference sample size;
- ancestry composition;
- genotype representation and allele counting;
- covariate adjustment;
- missing-data handling.

Do not treat these as incidental provenance. They determine whether the RSS
likelihood is compatible with the supplied statistics.

### 6. Handle residual variance according to the LD source

Estimating residual variance from RSS is recommended only when LD is genuinely
in sample and was constructed from the same effective genotype matrix used for
the marginal statistics. With external LD, residual variance should normally
remain fixed.

Current susieR also provides `R_finite` support for uncertainty caused by a
finite LD reference panel. When enabled, it reports diagnostics including
reference sample size, effective LD rank, effective-rank/reference-size ratio,
and per-variable penalties. The documentation states that an
effective-rank/reference-size ratio no greater than 0.2 indicates an adequate
reference panel under that diagnostic.

The 2026-07-28 official mismatch vignette additionally recommends considering
`R_mismatch = "eb"` for broader population/reference mismatch that remains
after finite-panel uncertainty is modeled. It must follow allele QC and should
be treated as a separately reviewed sensitivity when a primary workflow was
already fixed to `R_finite` alone.

Source:

- current `susie_rss()` reference:
  <https://stephenslab.github.io/susieR/reference/susie_rss.html>
- current external-LD mismatch vignette:
  <https://stephenslab.github.io/susieR/articles/rss_mismatch.html>

### 7. Run LD and summary-statistic mismatch diagnostics

At minimum run and record:

- `susieR::estimate_s_rss()`;
- `susieR::kriging_rss()`;
- observed versus expected z-score outliers;
- LD symmetry and unit diagonal;
- minimum and maximum eigenvalues;
- effective rank;
- variant ordering and allele concordance;
- genotype missingness, MAF, and monomorphic exclusions.

`estimate_s_rss()` estimates a regularization or inconsistency parameter. A
larger value indicates greater incompatibility between z-scores and the LD
matrix, but it is a diagnostic rather than a universal pass/fail statistic.
`kriging_rss()` helps localize variants with unexpected z-scores and is useful
for detecting allele flips or isolated mismatches.

Do not use these diagnostics to conceal a known allele or sample mismatch with
automatic correction. First correct harmonization, sample selection, and LD
construction.

Source:

- Stephens Lab RSS diagnostic vignette:
  <https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html>

### 8. Check optimization and credible-set quality

For every fit, report:

- convergence status;
- number of iterations;
- ELBO trajectory or final ELBO;
- maximum PIP and leading variants;
- number of non-null credible sets;
- credible-set membership and size;
- minimum, mean, and median credible-set correlation or purity;
- prior and residual variance estimates;
- warnings and numerical failures.

Slow or failed convergence can indicate unstable summary-statistic fitting,
LD mismatch, or a local optimum. `runsusie()` extends the run until convergence,
while direct `susie_rss()` requires the caller to manage iteration limits.

The `refine = TRUE` option can be used as a declared sensitivity for escaping
local optima. It should not substitute for fixing incompatible data.

Sources:

- current `susie_rss()` reference:
  <https://stephenslab.github.io/susieR/reference/susie_rss.html>
- coloc SuSiE vignette:
  <https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html>

### 9. Choose `L`, priors, and coverage explicitly

`L` is the maximum number of single-effect components, not the expected or
required number of causal signals. `L = 10` is a common default, but stability
should be checked when a region may contain many independent signals.

Record all effective parameters, including:

- `L`;
- `scaled_prior_variance`;
- whether prior variance is estimated;
- whether residual variance is estimated;
- credible-set coverage;
- minimum credible-set purity correlation;
- convergence tolerance and iteration limit;
- refinement setting;
- any finite-reference or mismatch model.

Do not reduce credible-set coverage simply to manufacture a signal. Absence of
a 95 percent credible set is a valid result.

### 10. Do not force numerical success by changing LD silently

Do not apply `nearPD`, arbitrary LD shrinkage, significance-based pruning, or
variant removal merely to make the algorithm converge. These operations alter
the statistical problem and may change PIPs and credible sets.

Exact or near-exact genotype ties should be retained when numerically possible.
If a computational representation must collapse tied variants, preserve a
complete reversible mapping back to every statistically indistinguishable
variant.

### 11. Validate coloc separately from SuSiE fitting

SuSiE fine-mapping output and `coloc.susie()` posterior probabilities are
different estimands and should not be numerically equal. Validation consists
of two separate checks:

1. Direct `susie_rss()` and `runsusie()` produce equivalent underlying SuSiE
   fits when inputs and arguments match.
2. Passing validated prefit SuSiE objects to `coloc.susie()` agrees with asking
   `coloc.susie()` to fit the same datasets internally with the same effective
   settings.

For colocalization, also report sensitivity to the cross-trait priors,
especially `p12`, and distinguish high H4 support from H3 support for separate
signals.

Sources:

- coloc SuSiE vignette:
  <https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html>
- coloc prior-sensitivity reference:
  <https://chr1swallace.github.io/coloc/reference/sensitivity.html>

## Recommendation for the targeted workflow

Use current `susieR` as the reference implementation and `coloc::runsusie()`
as the production coloc interface after equivalence is demonstrated. For the
MAPK3 pilot:

1. Fit direct individual-level `susie()` for each eQTL context using matched
   donors, genotypes, expression, and a rank-aware covariate model.
2. Fit direct `susie_rss()` using the stored dense tensorQTL beta and standard
   error statistics and matching covariate-adjusted in-sample eQTL LD.
3. Fit `runsusie()` from exactly the same z-scores, LD, sample size, variant
   order, and SuSiE parameters.
4. Compare PIPs, credible sets, log Bayes factors, ELBO, iterations,
   convergence, and ordering.
5. Run allele, LD, `estimate_s_rss()`, and `kriging_rss()` diagnostics before
   interpreting disagreement.
6. Fit the approved GWAS with ancestry-appropriate LD and explicitly account
   for its external finite-reference status.
7. Run `coloc.susie()` only after the trait-specific SuSiE fits pass review.
8. If all equivalence checks pass, use `runsusie()` for later targeted regions
   and parallelize independent CPU fits. Retain GPU tensorQTL as a possible
   separately validated scaling option, not as an untested replacement.

## Primary references

- Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to
  variable selection in regression, with application to genetic fine mapping.
  Journal of the Royal Statistical Society Series B. 2020.
  <https://doi.org/10.1111/rssb.12388>
- Zou Y, Carbonetto P, Wang G, Stephens M. Fine-mapping from summary data with
  the Sum of Single Effects model. PLOS Genetics. 2022.
  <https://doi.org/10.1371/journal.pgen.1010299>
- Wallace C. A more accurate method for colocalisation analysis allowing for
  multiple causal variants. PLOS Genetics. 2021.
  <https://doi.org/10.1371/journal.pgen.1009440>
- Taylor-Weiner A, Aguet F, Haradhvala NJ, et al. Scaling computational
  genomics to millions of individuals with GPUs. Genome Biology. 2019.
  <https://doi.org/10.1186/s13059-019-1836-7>
