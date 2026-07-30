# Workbook: 23andMe MDD and BD GWAS Integration and MBv Reanalysis

Audit date: 2026-07-30 EDT

Project repository:
`/home/gpertea/work/R/spatialDLPFC_mdd_bpd`

Project branch: `new-gwas`

GWAS integration repository:
`/dbdata/cdb/gwas-pgs-prs`

This document is a self-contained technical handoff for a reviewer with no
prior project context. It describes the delivered 23andMe data, the public GWAS
components, the reconstructed European meta-analyses, their limitations, and
the downstream MBv eQTL/GWAS and colocalization rerun.

## 1. Executive Summary

The approved 23andMe delivery did not contain final, already-integrated MDD or
BD European GWAS summary statistics. It contained the European 23andMe-only
association components used by the two papers, plus shared platform annotation
and QC tables.

For each disorder, the full European result was reconstructed from:

1. the existing public European GWAS that explicitly excluded 23andMe; and
2. the newly delivered European 23andMe-only association result.

The two component files were not concatenated. Shared variants were mapped to
exact alleles, normalized on GRCh37, and combined by standard-error
inverse-variance fixed-effect meta-analysis. This recalculated combined beta,
standard error, P value, odds ratio, effective sample size, and two-component
heterogeneity statistics. The result was then filtered by the paper-specific
effective-sample-size rule and lifted once to GRCh38.

The resulting files are usable for the MBv sensitivity analysis, PRS work, and
eQTL/GWAS colocalization, with these status labels:

- MDD: reconstructed 23andMe-inclusive European meta-analysis, release tag
  `full23andMe`.
- BD: reconstructed 23andMe-inclusive European meta-analysis, release tag
  `full23andMe_preDENTIST`.

The BD result must not be called the exact final paper GWAS. It uses a v7.0
association file with v7.2 annotations and lacks the paper's final
post-meta-analysis DENTIST exclusions. The MDD inputs are release-matched, but
the MDD file is still a reconstruction rather than a byte-identical official
final release.

Neither reconstructed BCF has a defensible combined imputation-quality score.
The public component supplies `INFO`; the 23andMe annotation supplies
`avg.rsqr` and `min.rsqr`. These are component-specific metrics and were not
averaged into a fabricated `SI`. Final `FORMAT/SI` values are missing.

The MBv `new-gwas` analysis then:

- rebuilt release-specific MDD and BD GWAS/genotype harmonization caches;
- reran strict and exploratory GWAS/eQTL/DEG overlap tables;
- reran `coloc::coloc.abf` and `coloc::sensitivity` for MDD and BD in all eight
  all-donor Seurat cell contexts;
- rebuilt aggregate coloc tables, figures, and final Excel workbooks; and
- compared the new results with a frozen no-23andMe archive.

The old results were preserved before current canonical output paths were
reused. Generated raw outputs are local ignored work products. Selected tables,
figures, documentation, and requested workbooks are committed on `new-gwas`.

This audit also identified a downstream metadata defect in the current coloc
runner. The integrated BCF query correctly reads effective N (`NE`), total
sample count (`NS`), and case count (`NC`), but the matched cache drops `NS`
and `coloc_case_fraction()` records `NC/NE` rather than the case fraction
`NC/NS`. The recorded `s` values are therefore wrong. This did not alter the
reported posterior probabilities in this run: `coloc` 5.2.3 uses the supplied
GWAS beta and variance directly and does not use `N` or `s` in that code path.
The metadata and cache schema should nevertheless be corrected before another
run or any switch to P-value/MAF-based coloc input.

## 2. Scope and Terminology

This audit covers only European-ancestry MDD and BD statistics. Non-European
BD files and the Hyde 2016 MDD delivery are out of scope.

Terminology used here:

- `public no23`: public European meta-analysis excluding 23andMe.
- `23andMe component`: delivered European association statistics generated
  within the 23andMe cohort.
- `integrated` or `full23`: the local public-no23 plus 23andMe reconstruction.
- `final paper GWAS`: the authors' final analysis after every paper-specific
  QC step. This term is not used for reconstructed BD.
- `significant variants`: variant rows at `P <= 5e-8`; these are not counts of
  independent signals or loci.
- `SI`: the GWAS-VCF field used for an imputation-accuracy statistic. It is
  missing in the reconstructed MDD and BD BCFs.

## 3. Study and Cohort Context

The relevant papers are:

- O'Connell et al., *Genomics yields biological and phenotypic insights into
  bipolar disorder*, Nature, DOI
  [10.1038/s41586-024-08468-9](https://doi.org/10.1038/s41586-024-08468-9).
- Adams et al., *Trans-ancestry genome-wide study of depression identifies 697
  associations implicating cell types and pharmacotherapies*, Cell, DOI
  [10.1016/j.cell.2024.12.002](https://doi.org/10.1016/j.cell.2024.12.002).

The cohort accounting used in this project is:

| Trait and component | Cases | Controls | N_eff/2 |
|---|---:|---:|---:|
| BD public EUR, no23 | 59,287 | 781,022 | 81,649 |
| BD delivered 23andMe EUR | 72,682 | 1,541,394 | 138,818 |
| BD paper full EUR | 131,969 | 2,322,416 | 220,467 |
| MDD public EUR, no23 | 412,305 | 1,588,397 | 576,327 |
| MDD delivered 23andMe EUR | 112,892 | 1,773,938 | 212,276 |
| MDD paper full EUR | 525,197 | 3,362,335 | 788,603 |

The delivered case/control counts came from the supplied `bipolar.html` and
`mdd.html` reports. They exactly equal the increments from each public-no23
cohort to its paper full-European cohort. This is strong evidence that the
delivered association files are the omitted 23andMe components, not separate
full meta-analyses.

The `N_eff/2` values above are paper-reported values. For an individual
case-control cohort, including each delivered 23andMe component:

```text
N_eff = 4 * N_case * N_control / (N_case + N_control)
```

The formula gives 138,818.2 for BD 23andMe and 212,275.0 for MDD 23andMe,
consistent with the reported increments after rounding. It must not be applied
to pooled case/control totals for a meta-analysis of cohorts with different
case/control ratios. Effective N is calculated within each contributing cohort
and then summed. This is why applying the formula to the displayed public or
full pooled counts gives different values. The integration script stores and
sums `N_eff`, not `N_eff/2`.

## 4. What the 23andMe Packet Contained

The delivery root is:

```text
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/
```

The relevant European association files are:

| Trait | Association file | Release |
|---|---|---:|
| BD | `Bipolar-Disorder-O_Connell-2025/OConnell_2025_bipolar_european-7.0/bipolar.dat.gz` | 7.0 |
| MDD | `MDD-Adams-2025/Adams_2025_mdd_european-7.2/mdd.dat.gz` | 7.2 |

Each association file has 57,611,376 data rows and this 16-column schema:

```text
all.data.id src pvalue effect stderr pass
im.num.0 dose.b.0 im.num.1 dose.b.1
AA.0 AB.0 BB.0 AA.1 AB.1 BB.1
```

Field interpretation from the supplied 23andMe documentation:

- `all.data.id`: internal integer key; it is not a genomic coordinate.
- `src`: `I` for imputed result or `G` for directly genotyped result.
- `pvalue`: component association-test P value.
- `effect`: log odds per copy of the B allele.
- `stderr`: standard error of `effect`.
- `pass`: 23andMe result-level QC status; only `Y` was retained.
- `.0` fields: controls.
- `.1` fields: cases.
- `im.num.*`: analyzed sample counts for imputed variants.
- `AA/AB/BB.*`: genotype counts used to reconstruct sample counts for
  genotyped variants.

The 23andMe A/B allele labels are alphabetical. The B allele is the
alphabetically higher allele. Coordinates and exact alleles are not present in
the association files and must be joined from the matching platform release.

The relevant European v7.2 annotation package contains:

```text
7.2-Annotations/v7.2_europe/all_snp_info.txt.gz
7.2-Annotations/v7.2_europe/gt_snp_stat.txt.gz
7.2-Annotations/v7.2_europe/im_snp_stat.txt.gz
7.2-Annotations/v7.2_europe/23andMe_GWAS_Results_v7.2.docx
7.2-Annotations/v7.2_europe/23andMe_Platform_Annotations_v7.2.docx
```

Validated annotation row counts are:

| File | Rows excluding header |
|---|---:|
| `all_snp_info.txt.gz` | 57,611,376 |
| `gt_snp_stat.txt.gz` | 1,481,446 |
| `im_snp_stat.txt.gz` | 57,525,634 |

`all_snp_info` maps `all.data.id` to GRCh37 scaffold, position, rsID/assay
name, alleles, ploidy, and strand. `im_snp_stat` includes `avg.rsqr`,
`min.rsqr`, `p.batch`, and `qc.mask`. `gt_snp_stat` includes genotype call
rate, frequencies, Hardy-Weinberg P value, and date-batch P value.

The supplied methods document states that imputed results with `rsq < 0.3`
were failed, as were results with strong platform-batch evidence. The
association `pass` flag therefore carries source QC decisions. When both
genotyped and imputed results passed, 23andMe generally reported the imputed
result.

Important limitation: the integration script consumed `all_snp_info` and the
association file. It did not join `im_snp_stat` or `gt_snp_stat` into the final
BCF. It relied on `pass=Y`, and it did not preserve `avg.rsqr`, `min.rsqr`, or
genotype call-rate fields as separate component provenance fields.

### 4.1 Packet integrity and inventory observations

Full-stream checks found:

| Check | BD component | MDD component |
|---|---:|---:|
| Total rows | 57,611,376 | 57,611,376 |
| Unexpected/nonsequential IDs | 0 | 0 |
| `pass=Y` rows with P | 23,844,053 | 23,795,959 |
| Raw `pass=Y`, `P <= 5e-8` rows | 4,855 | 12,568 |
| Passing imputed rows | 23,688,193 | 23,636,733 |
| Passing genotyped rows | 155,860 | 159,226 |

The vendor document describes a 14-column file, but the delivered header and
the document's own field table contain 16 columns. The delivered schema is
internally consistent; this appears to be a documentation counting error.

Original TAR files were retained. Extracted large annotation text files were
compressed individually as `.txt.gz`; gzip integrity tests passed. No original
archive was deleted.

### 4.2 Packet checksums used or recorded

| File | SHA-256 |
|---|---|
| BD association | `ef6647525f18a16cf585d4106179c3da233870385e08b989961c0e69dcb3878d` |
| MDD association | `29347df67ab0ff76dcff2c5ce1a0f3da90265e9d6d60fe15e2175f42c7bf92d6` |
| `all_snp_info.txt.gz` | `983428e8cced6200a2642f4a95629877ed0566d4fe7c823cecff57d2ea66937c` |
| `im_snp_stat.txt.gz` | `16c24cc54d9a118bfb1131cb40304e907c928a49f24c53fa741ab8f063056a0a` |
| `gt_snp_stat.txt.gz` | `1d572d89291e2e28b89b82d32e8130d25a23b86393c4000b1d2eccdaa876056c` |
| GWAS methods DOCX | `5d4fa461c03fa031bb2109e8fffd27b358aad686a7a04829164c2b9651210482` |
| Platform annotation DOCX | `90c7e0822e5cb59cb63f83a34d2ac00acac47ac513446613fdb5a24f81138a81` |

## 5. Public no-23andMe Components

The existing public European files were:

| Trait | Public source file | Rows | SHA-256 |
|---|---|---:|---|
| BD | `/home/gpertea/work/ref/GWAS/BD/bip2024_eur_no23andMe.gz` | 6,939,126 | `eccca2741f47416332f1a145aa3742367ec050a7ba9c6d7850ed3e5f06e06b13` |
| MDD | `/home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz` | 7,363,302 | `5d6fc5aee638e73457da703b75e0b4d4dacfb1b87578fd464b771209c0dfb22a` |

The public files already had GRCh37 coordinates, alleles, association
statistics, case/control counts, effective sample size, and an `INFO` field.
BD supplied odds ratios, which `bcftools +munge` converted to log effects. MDD
supplied beta directly.

The public files are available through:

- BD Figshare: [10.6084/m9.figshare.27216117](https://doi.org/10.6084/m9.figshare.27216117).
- MDD Figshare: [10.6084/m9.figshare.27061255](https://doi.org/10.6084/m9.figshare.27061255).

The public files exclude the delivered 23andMe cohorts. This exclusion and the
exact case/control increments support, but do not participant-level prove, the
no-overlap assumption required by the fixed-effect integration.

## 6. Why Concatenation Would Be Wrong

Each source row is an association estimate from a different participant set.
Concatenating rows would create duplicate variant records and would not combine
evidence. Adding P values or selecting the smaller P value would also be
invalid.

For public component `p` and 23andMe component `m`, the integration uses:

```text
w_p = 1 / SE_p^2
w_m = 1 / SE_m^2

BETA_full = (w_p * BETA_p + w_m * BETA_m) / (w_p + w_m)
SE_full   = sqrt(1 / (w_p + w_m))
Z_full    = BETA_full / SE_full
P_full    = 2 * Phi(-abs(Z_full))
OR_full   = exp(BETA_full)
```

This is standard standard-error inverse-variance fixed-effect meta-analysis.
It is appropriate if:

- both components estimate the same log-odds effect;
- effect alleles are aligned;
- standard errors are valid;
- component participant sets are independent; and
- a common underlying fixed effect is the intended analysis model.

A prior public meta-analysis can be combined with an independent aggregate
component because fixed-effect inverse-variance meta-analysis is associative.
Individual-level genotypes are not required.

The method is implemented by `bcftools +metal`, following the METAL
inverse-variance approach, rather than by custom arithmetic. Primary method
reference: [Willer, Li, and Abecasis 2010](https://doi.org/10.1093/bioinformatics/btq340).
GWAS meta-analysis QC reference:
[Winkler et al. 2014](https://doi.org/10.1038/nprot.2014.071).

## 7. Integration Software and Repository Provenance

Reusable integration code is in the `gwas-pgs-prs` repository:

```text
/dbdata/cdb/gwas-pgs-prs/mbv-prs/scripts/integrate_23andme_gwas.sh
/dbdata/cdb/gwas-pgs-prs/mbv-prs/scripts/prepare_23andme_gwas.awk
/dbdata/cdb/gwas-pgs-prs/mbv-prs/scripts/gwas_meta_colheaders.tsv
/dbdata/cdb/gwas-pgs-prs/mbv-prs/scripts/gwas_meta_extra_headers.txt
/dbdata/cdb/gwas-pgs-prs/mbv-prs/tests/test_gwas_meta_integration.sh
```

Method documentation is in:

```text
/dbdata/cdb/gwas-pgs-prs/mbv-prs/GWAS_23andMe_integration.md
/dbdata/cdb/gwas-pgs-prs/mbv-prs/BD-DENTIST-issue-v7.0-vs-v7.2.md
```

Current integration-repository commit at audit time:

```text
5d3db2c18894c901d6851f30dfe9d63de6a535cd
```

The integration workflow was introduced in commit `1e5d7b8` and is present in
the current `master` commit above.

The production run reports did not record executable versions or hashes.
Current srv16 versions and hashes were therefore captured on 2026-07-30 as
audit context, not asserted as immutable historical run metadata:

| Item | Current audit value |
|---|---|
| `bcftools` | 1.23.1, htslib 1.23.1 |
| `/opt/sw/bin/bcftools` SHA-256 | `318d2784eeb97c56e7df325fc82e295653f0f5454ebfa8462288c1b2802684ba` |
| `munge.so` SHA-256 | `8bb3088a32e9d0781dbed61c95069b1bd81f288c70cefaf139232b716715eef3` |
| `metal.so` SHA-256 | `bf15bcfc19e218f41c59b191e11a909b11ed486fff4c9b23fe824b8edcf1ee69` |
| `liftover.so` SHA-256 | `13d90239a2ac95fe82fd6c51a66fd9df13d9f9789345d16322ac2f4ecfa7153b` |
| Plugin directory | `/opt/sw/bcf-plugins` |

## 8. Reference Genome and Liftover Inputs

The production script used:

| Purpose | Path | SHA-256 |
|---|---|---|
| GRCh37 normalization | `/dbdata/cdb/ref/GRCh37/human_g1k_v37.fasta` | `2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f` |
| GRCh38 target | `/dbdata/cdb/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7` |
| GRCh37/hg19 to GRCh38 chain | `/dbdata/cdb/ref/hg19ToHg38.over.chain.gz` | `5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d` |

The user also maintains reference resources under `~/work/ref/`, including
GRCh38 and 1000 Genomes resources. Reviewers should inspect that location first
for future LD or genome-reference needs. The completed integration itself used
the exact paths and checksums above.

## 9. Detailed Integration Pipeline

### Step 1: Fail-fast input and environment validation

The shell runner requires:

- trait exactly `BD` or `MDD`;
- an output prefix;
- an explicit `--confirm-no-sample-overlap` assertion;
- readable, nonempty source/reference/helper files;
- positive thread count;
- `awk`, `bcftools`, `bgzip`, `gzip`, `paste`, `sha256sum`, and `tabix`;
- `bcftools` plugins `munge`, `metal`, and `liftover`.

Existing outputs cause failure unless `--resume` is supplied. `--resume` only
reuses nonempty stage files; it is not a force-overwrite mode. This protects
large intermediates from accidental replacement.

Justification: source mix-ups or partial reruns are more damaging than a hard
failure. An explicit no-overlap assertion makes the key statistical assumption
visible at invocation time.

### Step 2: Stream annotation and association rows together

The script uses `paste` on decompressed `all_snp_info` and association streams,
then passes 34 combined fields to `prepare_23andme_gwas.awk`.

The AWK parser verifies:

- expected annotation and association headers;
- exactly 34 pasted fields;
- equality of the annotation and association `all.data.id` values;
- strict sequential IDs matching the data-row number.

Any ID mismatch causes nonzero exit after reporting examples. Production runs
had zero ID errors across 57,611,376 rows for both traits.

Justification: the association file has no coordinates or explicit alleles. A
one-row offset would silently assign every statistic to the wrong variant.

### Step 3: Apply 23andMe result and variant filters

The parser retains rows satisfying all of the following:

- `pass=Y`;
- chromosome 1 through 22;
- autosomal ploidy code;
- forward strand;
- exactly two one-base A/C/G/T alleles;
- numeric beta, positive SE, and valid P in `(0,1]`;
- positive case and control counts;
- recognized source `I` or `G`.

Imputed sample counts come from `im.num.0` and `im.num.1`. Genotyped counts are
the sums of AA, AB, and BB counts by phenotype group.

The effect allele is emitted as alphabetical B, and the other allele as A.
Per-row effective N is calculated from the reconstructed case/control counts.

Production parser counts were:

| Filter outcome | BD | MDD |
|---|---:|---:|
| Total input rows | 57,611,376 | 57,611,376 |
| Kept autosomal SNVs | 21,137,709 | 21,098,300 |
| Failed 23andMe QC | 33,767,323 | 33,815,417 |
| Symbolic/non-SNV | 2,056,346 | 2,051,820 |
| Nonautosomal | 649,998 | 645,839 |
| Invalid statistics | 0 | 0 |
| Invalid sample size | 0 | 0 |
| Invalid source | 0 | 0 |
| ID errors | 0 | 0 |

Justification: exact sequence alleles are required for safe normalization,
liftover, target-genotype matching, and beta orientation. Symbolic `D/I`
alleles do not identify an insertion/deletion sequence and were excluded.

### Step 4: Convert both components to normalized GRCh37 BCF

`bcftools +munge` uses `gwas_meta_colheaders.tsv` to recognize alternate source
column names. Important mappings include:

```text
SNP/ID -> SNP
CHR/CHROM/#CHROM -> CHR
BP/POS -> BP
A1/EA -> A1
A2/NEA -> A2
OR -> OR
BETA -> BETA
SE -> SE
P/PVAL -> P
INFO/IMPINFO -> INFO
NCAS/Nca -> N_CAS
NCON/Nco -> N_CON
NEFF -> NEFF
NEFFDIV2/Neff_half -> NEFFDIV2
```

The public and prepared 23andMe tables are independently converted using the
same synthetic trait sample name. BD public OR is converted to log effect;
MDD and 23andMe beta stay on the log-odds scale.

Both streams then pass through:

```text
bcftools view -f PASS,. -m2 -M2 -v snps
bcftools norm -f GRCh37 -c e -d exact
bcftools sort
```

This retains biallelic SNVs, validates reference alleles, removes exact
duplicates, and produces sorted indexed BCFs.

Normalization removed 41 exact duplicates from the BD 23andMe component and
49 from MDD. No GRCh37 reference mismatches were reported in either public
component.

Justification: rsIDs are neither unique nor stable enough for allele-sensitive
meta-analysis. GRCh37 position plus normalized REF/ALT defines the join.

### Step 5: Run fixed-effect inverse-variance meta-analysis

The command is:

```text
bcftools +metal --het --esd PUBLIC.bcf 23ANDME.bcf
```

`--esd` selects standard-error/effect-size meta-analysis. `--het` emits
two-component Cochran heterogeneity fields. Alleles are aligned by the BCF
representation; effect direction is flipped where required.

The meta-analysis outputs:

- `ES`: combined log-odds beta relative to ALT;
- `SE`: combined standard error;
- `LP`: combined `-log10(P)`;
- `NS`: summed analyzed sample count;
- `NC`: summed case count;
- `NE`: summed effective sample size;
- `I2`: Cochran I-squared across the two aggregates;
- `CQ`: Cochran Q `-log10(P)`;
- `ED`: effect directions across the two aggregates.

The delivered component P value is not averaged with the public P value. The
combined P follows from combined beta and SE.

Justification: both papers used fixed-effect inverse-variance meta-analysis,
and the component effect and SE fields are available. This is preferable to a
sample-size-weighted Z method because interpretable effect estimates and
standard errors are available.

### Step 6: Do not apply an additional genomic-control multiplier

The 23andMe methods document says raw association files are not adjusted for
genomic inflation. The supplied HTML displays lambda values of 1.122 for BD and
1.222 for MDD. That does not establish that the raw SE should be multiplied by
`sqrt(lambda)` before reconstructing the paper result.

This was checked against published full-European sentinel effects:

| Trait | Sentinels | Mean absolute beta error, raw SE | Error after extra lambda correction |
|---|---:|---:|---:|
| BD | 258 | `3.06e-5` | `3.52e-4` |
| MDD | 615 | `3.96e-5` | `3.02e-4` |

The raw beta/SE integration matched published values materially better. No
extra genomic-control correction was applied.

Justification: applying an unsupported second correction would downweight the
23andMe component and move reconstructed effects away from published results.

### Step 7: Apply paper-specific effective-N coverage filtering

The script finds maximum `NE` after meta-analysis and retains:

- BD: `NE >= 0.75 * max(NE)`;
- MDD: `NE >= 0.80 * max(NE)`.

Production thresholds were:

| Trait | Observed maximum NE | Expected paper-scale maximum | Ratio | Fraction | Minimum NE |
|---|---:|---:|---:|---:|---:|
| BD | 440,999 | 440,934 | 1.00015 | 0.75 | 330,749.25 |
| MDD | 1,577,200 | 1,577,206 | 1.00000 | 0.80 | 1,261,760 |

The near-exact maximum-NE agreement is a strong sample-accounting check.

This threshold also explains why the final files have about 6.4 to 6.7 million
variants rather than all 21 million passing 23andMe autosomal SNVs. Variants
available in only one component generally do not reach the required fraction
of full effective N.

Justification: per-variant missingness and cohort coverage vary. The threshold
implements each paper's coverage rule and avoids treating low-coverage variants
as if they represented the full analysis.

### Step 8: Optional post-meta exclusion list

The runner supports `--exclude-list` after effective-N filtering. No exclusion
list was supplied for either production run.

For MDD, no additional unavailable paper-specific exclusion list was identified
in this audit. For BD, absence of the final DENTIST exclusions is a known
material limitation described separately below.

### Step 9: Lift the combined result once to GRCh38

Only the final filtered GRCh37 meta-analysis is lifted:

```text
bcftools +liftover -s GRCh37.fa -f GRCh38.fa -c hg19ToHg38.chain
bcftools annotate -h gwas_meta_extra_headers.txt
bcftools sort
```

The integration does not independently lift each component before combining.
This reduces opportunities for build-specific join loss.

Liftover audit counts:

| Trait | Input lines | Allele-swapped | Reference added | Rejected | Final records |
|---|---:|---:|---:|---:|---:|
| BD | 6,394,928 | 8,228 | 5 | 140 | 6,394,788 |
| MDD | 6,656,330 | 7,799 | 3 | 108 | 6,656,222 |

`ES` and `ED` are handled with allele-flip rules during liftover; additive
fields are retained under the plugin's aggregation rules.

Justification: effect orientation must follow the target-build ALT allele. A
coordinate-only liftover without REF/ALT and effect handling would be unsafe.

### Step 10: Represent unavailable combined imputation quality honestly

The final header adds:

```text
FORMAT/SI: Combined imputation quality is unavailable for the reconstructed
public plus 23andMe meta-analysis
```

Every integrated `SI` value is missing.

Why:

- public `INFO` describes public-component imputation quality;
- 23andMe `avg.rsqr` and `min.rsqr` describe its platform/imputation component;
- `bcftools +metal` does not emit a combined imputation score; and
- no justified arithmetic rule was identified that would reproduce a
  paper-equivalent combined quality statistic.

This does not mean every variant has poor imputation. It means the reconstructed
file cannot support one post-integration `SI` threshold.

Reviewer concern: the integration script also did not retain the two component
metrics in separate fields. A future revision should preserve, for example,
`PUBLIC_INFO`, `23ME_AVG_RSQ`, and `23ME_MIN_RSQ` or a sidecar table. Those
fields must not be relabeled as combined `SI`.

### Step 11: Export reusable products

For each prefix, the runner writes:

| Product | Purpose |
|---|---|
| `PREFIX.23andme.grch37.ssf.tsv.gz` | Prepared 23andMe component |
| `PREFIX.public-no23.grch37.bcf` | Normalized public component |
| `PREFIX.23andme.grch37.bcf` | Normalized 23andMe component |
| `PREFIX.meta-unfiltered.grch37.bcf` | Unfiltered two-component meta-analysis |
| `PREFIX.full.grch37.bcf` | Effective-N-filtered GRCh37 meta-analysis |
| `PREFIX.full.grch38.bcf` | Dense GRCh38 analysis BCF |
| `PREFIX.full.grch38.meta.tsv.gz` | Canonical GRCh38 summary-statistics registry |
| `PREFIX.full.grch38.prsice.tsv.gz` | PRSice base table |
| `PREFIX.integration-report.md` | Input hashes, parameters, counts, caveats |

The canonical table fields are:

```text
#CHROM POS RSID VARIANT REF ALT EA NEA BETA SE P OR NS NC NE
```

`VARIANT` is `CHROM:POS:REF:ALT`; EA is ALT. The PRSice table uses:

```text
CHR BP SNP A1 A2 BETA P
```

with canonical SNP ID `CHROM:POS:REF:ALT`, A1=ALT, and A2=REF.

## 10. Production Commands

MDD used matching v7.2 association and annotation releases:

```bash
cd /dbdata/cdb/gwas-pgs-prs/mbv-prs

scripts/integrate_23andme_gwas.sh \
  --trait MDD \
  --out-prefix /home/gpertea/work/ref/GWAS/MDD/full_eur_integration/pgc-mdd2025_eur_v3-49-24-11 \
  --confirm-no-sample-overlap \
  --threads 8
```

BD required explicit acknowledgement of the v7.0/v7.2 mismatch:

```bash
cd /dbdata/cdb/gwas-pgs-prs/mbv-prs

scripts/integrate_23andme_gwas.sh \
  --trait BD \
  --out-prefix /home/gpertea/work/ref/GWAS/BD/full_eur_integration/bip2024_eur \
  --confirm-no-sample-overlap \
  --allow-bd-v7.2-annotations \
  --threads 8
```

The completed logs are:

```text
/home/gpertea/work/ref/GWAS/BD/full_eur_integration/integration.run.log
/home/gpertea/work/ref/GWAS/MDD/full_eur_integration/integration.run.log
```

Run timestamps:

- BD: 2026-07-26 17:27:33 through 17:43:36 EDT.
- MDD: 2026-07-26 17:27:41 through 17:40:42 EDT.

## 11. Integrated Production Artifacts

Stable installed products omit `_no23andMe` from their names:

| Trait | Product | Size | SHA-256 |
|---|---|---:|---|
| BD | `BD/bip2024_eur.hg38.bcf` | 213 MB | `1d502351659d81aa503da64779b7cc9e80c6819bca9e2ab5d5ac6efc3bb531be` |
| BD | `BD/bip2024_eur.hg38.meta.tsv.gz` | 223 MB | `6c0f7b2075407085ba490cb6dd0d64bb03c3bc64ca1dc8df149919dcef29cb27` |
| BD | `BD/bip2024_eur.hg38.prsice.tsv.gz` | 112 MB | `fb32072d7bc035c0d3e8dcc6232cefdf5a68a97ffbe2c7cbdedf5cadf8cf377c` |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf` | 222 MB | `a8b30df37b032097920ded06697aa60851f0078e681f1a595b6d5f3ec70347f3` |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.meta.tsv.gz` | 233 MB | `53cd51f010bc1146d9635eb9170a70e48d3ecc8018d2a58b16a9a3dc0b3f3680` |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.prsice.tsv.gz` | 118 MB | `d4d4a5d457a11d794aa3f03c3043afa112edd368cdf46b0ab9b22124287b5b43` |

The stable paths under `/home/gpertea/work/ref/GWAS` resolve to shared storage
under `/dbdata/cdb/ref/GWAS`. The spatialDLPFC project accesses the same files
through `processed-data/ref/GWAS`.

Final BCF summary:

| Trait | Synthetic sample | Records | `P <= 5e-8` variants |
|---|---|---:|---:|
| BD | `BD_2024_FULL_EUR` | 6,394,788 | 10,753 |
| MDD | `MDD_2025_FULL_EUR` | 6,656,222 | 36,787 |

The final FORMAT fields are:

```text
NS:NC:ES:SE:LP:NE:I2:CQ:ED
```

`SI` is defined in the header but values query as missing.

## 12. Trait-Specific Caveats

### 12.1 BD association v7.0 with annotation v7.2

The BD TAR contains the association file and HTML reports, but no European v7.0
`all_snp_info`, `gt_snp_stat`, or `im_snp_stat`. The only delivered European
annotation package is v7.2.

Both BD v7.0 association and v7.2 annotation contain the same sequential
57,611,376-ID universe. Published-sentinel agreement supports the mapping at
checked loci. However, the supplied documents do not guarantee stable
`all.data.id`, coordinates, or alleles across releases.

Consequence: genome-wide BD mapping remains provisional. A reviewer should
seek either:

1. the correct European v7.0 annotation package;
2. written 23andMe confirmation of v7.0-to-v7.2 ID, coordinate, and allele
   stability; or
3. a complete validated cross-release map.

### 12.2 BD post-meta DENTIST not reproduced

The public no23 BD file was already DENTIST-filtered. The paper then applied
DENTIST after the full meta-analysis using ancestry-matched HRC LD. The local
reconstruction created new combined z scores, but no final-paper DENTIST
exclusion list, exact software version, parameters, HRC files, or matching rules
were delivered.

Filtering the public component before integration does not substitute for
post-meta QC of newly combined statistics.

Consequence: label the result `full23andMe_preDENTIST`. The runner supports a
paper-provided GRCh37 exclusion list via `--exclude-list`, but none was
available.

References:

- [BD paper methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12163093/).
- [DENTIST software](https://github.com/Yves-CHEN/DENTIST).

### 12.3 MDD status

MDD uses matching v7.2 association and annotation releases, and the 80%
effective-N rule is documented in the public file metadata and paper methods.
Maximum effective N and 615 sentinel effects agree closely with published
statistics.

Consequence: MDD has substantially stronger reconstruction provenance than BD.
It is still locally reconstructed rather than an official final file supplied
as one artifact.

## 13. Cross-Trait Caveats and Assumptions

### 13.1 Participant independence

Fixed-effect standard errors are valid only if the public and 23andMe
components are independent or covariance is modeled. The script requires the
operator to assert no overlap. Evidence is the public no23 definition and exact
cohort-count increments. No participant-level crosswalk was available.

### 13.2 Fixed-effect model

The model assumes one common underlying effect per variant. `I2`, `CQ`, and
`ED` can flag disagreement between the two aggregates, but no heterogeneity
filter was applied. The emitted heterogeneity is between two aggregate
components, not among all original cohorts.

### 13.3 Imputation quality and downstream SI filtering

23andMe did provide component QC data, but the integration did not preserve it
as separate final fields. The public-no23 spatial workflow previously applied
`SI >= 0.8`. Integrated MDD/BD has missing SI and therefore receives no
post-integration SI filter.

The integration script also did not explicitly prefilter public rows at
`INFO >= 0.8`; it retained source PASS/unfiltered rows and later applied the
effective-N threshold. Consequently, downstream differences are not a pure
counterfactual that changes only participant count. They also reflect the
integrated file's row-level QC representation and absence of the previous
post-query SI threshold.

### 13.4 Symbolic indels

All 23andMe D/I variants were excluded because sequence alleles were unknown.
The public files contain only two BD and one MDD indels, so this does not remove
a substantial shared variant class. It nevertheless prevents claiming exact
all-variant paper reproduction.

### 13.5 Genomic inflation

No extra lambda correction was applied because sentinel comparisons strongly
favored raw delivered beta/SE. If an official full file becomes available,
this decision should be checked again directly.

### 13.6 P values and precision

Final P values were recalculated from combined beta and SE. Exported tables
derive P from BCF `LP`; underflow is capped at `1e-300`. Text exports use finite
decimal precision, while BCF remains the preferred dense analysis source.

### 13.7 European-only scope

No non-European component was integrated. Conclusions here concern European
summary statistics and the MBv target genotype panel.

## 14. Integration Validation

### 14.1 Unit fixture

`tests/test_gwas_meta_integration.sh` creates small annotation, association,
and public fixtures. It verifies:

- parser/header behavior;
- effect-allele orientation;
- exact and swapped allele handling;
- expected combined beta;
- expected combined SE; and
- summed effective N.

It was rerun successfully on 2026-07-30:

```text
PASS: parser, allele orientation, beta, SE, and effective N
```

### 14.2 Production checks

Verified production evidence includes:

- zero annotation/association ID alignment errors;
- zero invalid statistic/sample/source rows among otherwise retained records;
- maximum NE within 0.015% of expected paper scale for BD and effectively
  identical for MDD;
- superior sentinel agreement without extra lambda correction;
- readable BCF indexes reporting exactly 6,394,788 BD and 6,656,222 MDD rows;
- explicit input SHA-256 values in each integration report; and
- preserved normalized component and unfiltered meta-analysis BCFs for audit.

## 15. Activation in the spatialDLPFC MBv Project

The downstream project is:

```text
/home/gpertea/work/R/spatialDLPFC_mdd_bpd
```

The integration reanalysis is on branch `new-gwas`, based on `devel` commit:

```text
4bd471a64dca1cf663ecdc18e0d6050ac3a07f1b
```

Relevant branch commits are:

| Commit | Purpose |
|---|---|
| `56814982` | Switch GWAS inputs, add release-aware caches/tests, rerun code/docs/tables |
| `4991e5aa` | Commit five explicitly requested regenerated PDF/XLSX outputs |
| `ad11d80a` | Document integrated-GWAS coloc output locations and naming policy |

The impact report is committed at:

```text
code/11_eQTL_coloc/GWAS-23andMe-MDD-BD-impact-on-tables.md
```

### 15.1 GWAS selection and provenance in `utils.R`

`code/11_eQTL_coloc/utils.R` now selects:

```text
BD  -> BD/bip2024_eur.hg38.bcf
MDD -> MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf
SCZD -> unchanged public SCZD BCF
```

Release tags are:

```text
BD = full23andMe_preDENTIST
MDD = full23andMe
SCZD = public
```

SI rules are:

```text
BD = none
MDD = none
SCZD = SI >= 0.8
```

The code calculates and records GWAS BCF SHA-256 values. Final Excel workbooks
include a `gwas_provenance` sheet containing release, BCF filename/hash, SI
rule, thresholds, and matched variant counts.

### 15.2 Release-specific cache names

Cache names include release and SI policy. Examples:

```text
GWAS-BD_full23andMe_preDENTIST_flt_p1e-5_SInone_...
GWAS-MDD_full23andMe_flt_p1e-5_SInone_...
gwas_astro_full23andMe_preDENTIST_nominal-variants_SInone.tsv.gz
gwas_astro_full23andMe_nominal-variants_SInone.tsv.gz
```

Justification: old no23 cache names must not be silently reused after switching
GWAS releases. Cache names are provenance controls, not cosmetic labels.

### 15.3 MBv target genotype harmonization

The fixed all-donor target is:

```text
processed-data/00_genotypes/plink2/merged_maf05
```

Audit values:

| Item | Value |
|---|---:|
| Donors in PSAM | 119 |
| Variants in PVAR | 6,056,513 |
| PSAM SHA-256 | `42293f81f12ce7d078748d0f7e901b12795b7482b3894d13e06e618fc33c8cdb` |
| PVAR SHA-256 | `fd2c2ea9c92b7d3923567a2fbfa3e38e87851693591982c4f9b8bca70b3bab76` |

Sparse significance/overlap caches and dense coloc caches both match GWAS to
the PVAR using GRCh38 chromosome, position, REF, and ALT.

- exact REF/ALT matches retain beta;
- swapped REF/ALT matches negate beta;
- exact matches are preferred if duplicate target IDs arise;
- rsID alone is never the harmonization key;
- the resulting `variant_id` is the PVAR canonical ID.

The R test fixture
`code/11_eQTL_coloc/tests/testthat/test-gwas-harmonization.R` verifies cache
names, BCF path resolution, GWAS query normalization, cached and direct query
loading, and exact/swapped/unmatched target handling. It passed on 2026-07-30.

## 16. Downstream Workload Executed

### 16.1 Freeze the no23 state

Before rebuilding current canonical outputs, the previous no23 MDD/BD results
were preserved at:

```text
processed-data/11_eQTL_coloc/seurat/archive/no23andMe_2026-07-28/
```

The archive contains:

- old MDD/BD per-context coloc outputs;
- old dense GWAS caches;
- old large `.qs2` coloc result objects;
- old aggregate coloc tables;
- old 03-series eQTL/GWAS tables;
- old final Excel workbooks;
- old plots;
- GWAS and gene-list checksums; and
- a 143-row pre-move path/size/mtime manifest.

The exact archive shell transcript was not retained. The manifest provides
path, size, and modification-time evidence, but not a cryptographic checksum
for every archived generated file. This is a provenance gap.

### 16.2 Rerun 03-series eQTL/GWAS/DEG overlaps

The following were executed using the integrated GWAS configuration:

```text
03_eqtl_explore.Rmd
03a_nominal_eQTLs.Rmd
03b_eQTL_boxplots.Rmd
03c_GWAS_relaxed_eQTLs.R
```

The strict exact-variant threshold remained `P <= 5e-8`. The exploratory MDD
and BD threshold remained `P < 1e-5`. Curated paper gene lists were held fixed.
Only exact canonical target variants were used for variant overlap.

The run logs and rendered audit notebooks are under:

```text
processed-data/11_eQTL_coloc/seurat/run_logs/23andMe_2026-07-28/
```

### 16.3 Rebuild dense coloc GWAS caches

For each disorder and cell context, `04_run_coloc.R`:

1. read full tensorQTL nominal cis parquet data;
2. converted nominal target variant IDs to one-base query regions using the
   PVAR as coordinate truth;
3. queried the full integrated BCF at all nominal eQTL positions;
4. applied no GWAS P-value filter;
5. applied no SI filter for MDD/BD because SI is unavailable;
6. harmonized exact and swapped alleles to PVAR;
7. flipped beta for swapped alleles; and
8. wrote release-specific dense caches.

No GWAS-significance prefilter is appropriate for coloc. Coloc needs dense
regional association statistics to compare causal hypotheses and LD-correlated
signals.

### 16.4 Rerun `coloc::coloc.abf`

`coloc::coloc.abf` was rerun for MDD and BD, for all eight all-donor Seurat
contexts:

```text
Astro, Inhb, L2.3, L4, L5, L6, Micro.Vasc, Oligo
```

SCZD was not rerun; its existing results were retained for aggregate tables.

Per-gene input construction as actually executed:

- eQTL dataset type: quantitative;
- eQTL beta/variance/P/MAF from tensorQTL nominal results;
- eQTL N: 119 from the preparation manifest;
- GWAS dataset type: case-control;
- GWAS beta/variance/P and variant-specific effective N from the integrated
  BCF;
- GWAS `N` supplied to coloc: median valid variant-specific `NE` per locus;
- GWAS `s` supplied to coloc: median `NC/NE`, due to the implementation defect
  described below, rather than the required `NC/NS` case fraction;
- minimum overlapping SNPs: 10;
- minimum eQTL evidence: at least one `abs(beta/SE) >= 2`;
- cis span: tensorQTL nominal 1 Mb window;
- priors: `p1=1e-4`, `p2=1e-4`, `p12=1e-5`.

The implementation calls `coloc::check_dataset` before `coloc::coloc.abf`.
Errors are captured per gene rather than terminating all contexts.

Production execution used 12 `BiocParallel::MulticoreParam` workers per
disorder run. BLAS, OpenMP, Arrow, and data.table worker threads were capped at
one to avoid nested oversubscription.

Run summary:

| Trait | Contexts | Gene inputs per context | Gene errors | Candidate rows | Sum of per-context elapsed time |
|---|---:|---:|---:|---:|---:|
| BD | 8 | 8,102 to 14,834 | 0 | 679 | 2.61 hours |
| MDD | 8 | 8,102 to 14,838 | 0 | 1,280 | 2.57 hours |

The two disorder runs overlapped in wall-clock time. Logs span approximately
23:15 EDT on 2026-07-28 to 01:52 EDT on 2026-07-29.

#### 16.4.1 Coloc case-fraction audit finding

The BCF query emits fields in this order:

```text
NE NS NC
```

`normalize_coloc_gwas_query_table()` names those values `N`, `ns`, and `ncas`.
`match_coloc_gwas_geno()` then retains `N` and `ncas` but drops `ns`.
Finally, `coloc_case_fraction()` calculates `ncas/N`. Thus the executed value
is `NC/NE`, despite the function comment calling it a case fraction.

Observed values were:

| Trait | Recorded median `NC/NE` | Correct full-BCF median `NC/NS` |
|---|---:|---:|
| BD | 0.2992494 | 0.0537664 |
| MDD | 0.3329933 | 0.1350979 |

The numerical coloc results do not need to be withdrawn solely because of this
defect. The installed `coloc` version is 5.2.3. Inspection of its
`process.dataset()` implementation confirms that when `beta` and `varbeta`
are supplied, as they were here, the case-control approximate Bayes factor is
calculated from beta, variance, and the fixed case-control prior SD. `N` and
`s` are not used. `check_dataset()` only checked that the supplied `s` was
between zero and one. Therefore correcting `s` metadata without changing any
other input should reproduce the current PP values exactly.

Required code follow-up:

1. retain `ns` through `match_coloc_gwas_geno()` and the dense cache;
2. calculate case fraction as `ncas/ns`;
3. record `NE`, `NS`, and `NC` with unambiguous names rather than aliasing
   `NE` as `N`;
4. add a test that distinguishes `NC/NS` from `NC/NE`; and
5. decide explicitly whether any future P-value/MAF coloc mode should use total
   `NS` or an effective-N approximation for `N`.

Until fixed, current run metadata must not be cited as evidence of the true
case fraction. The PP4 comparison remains valid for the current beta/SE-based
implementation, but this exception is version- and code-path-specific.

### 16.5 Coloc categorization and sensitivity gate

Raw categories are assigned from posterior probabilities:

```text
strong_coloc:   PP4 > 0.8
likely_coloc:   PP3+PP4 > 0.8 and PP4/(PP3+PP4) >= 0.9
distinct_causal: PP3+PP4 > 0.8 and PP4/(PP3+PP4) < 0.5
follow_up:      PP3+PP4 > 0.8 or PP4/(PP3+PP4) > 0.8
```

`coloc::sensitivity` evaluates rule `H4 > 0.8` on a 100-point prior grid. The
plausible `p12` interval is `5e-6` through `5e-5`.

A final gated strong call requires:

- raw category `strong_coloc`;
- no sensitivity error;
- at least one plausible grid point; and
- at least 50% of plausible points passing the rule.

Justification: PP4 at one chosen prior can be unstable. The gate requires a
strong call to persist over a plausible shared-causal-prior range.

### 16.6 Lean coloc output decision

The new runs did not use `--write-full-qs2`. Full in-memory `coloc.abf` objects
were therefore not persisted for the new MDD/BD runs.

Each context instead retains:

```text
coloc_<dataset>.flat.tsv.gz
coloc_<dataset>.sensitivity.tsv.gz
coloc_<dataset>.errors.tsv.gz
coloc_<dataset>.runmeta.tsv.gz
coloc_<dataset>.complete
```

The `.flat.tsv.gz` files contain candidate-level posterior summaries, not all
per-SNP posterior rows. The sensitivity files contain per-gene summarized
sensitivity-grid results. The old no23 full `.qs2` objects remain in the dated
archive.

Consequence: current candidate/final reporting is reproducible from persisted
lean outputs, but retrieving arbitrary new SNP-level posterior details from the
new full result objects would require rerunning with `--write-full-qs2` or
adding a deliberate SNP-level export.

### 16.7 Aggregate and final-table reruns

After the two `04` runs:

- `05_coloc_explore.Rmd` rebuilt candidate, gated, sensitivity, QC, summary,
  workbook, and plot outputs.
- `07_compare_23andme_impact.R` compared current results to the frozen no23
  archive and wrote the impact report and machine-readable deltas.
- `06_eqtl_coloc_final_tables.Rmd` rebuilt manuscript-facing eQTL and coloc
  workbooks with GWAS provenance sheets.

All 16 current MDD/BD context completion markers exist. Four principal XLSX
containers were retested with `unzip -t` on 2026-07-30.

## 17. Downstream Output Locations

### Current per-context and cache outputs

```text
processed-data/11_eQTL_coloc/seurat/coloc/BD/
processed-data/11_eQTL_coloc/seurat/coloc/MDD/
```

### Current aggregate coloc outputs

```text
processed-data/11_eQTL_coloc/seurat/coloc/tables/coloc_abf_results.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/tables/coloc_abf_results_gated.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/tables/coloc_abf_sensitivity.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/tables/coloc_abf_qc.tsv.gz
processed-data/11_eQTL_coloc/seurat/coloc/tables/coloc_abf_results.xlsx
```

### Current final workbooks

```text
processed-data/11_eQTL_coloc/seurat/final/eqtl_results.xlsx
processed-data/11_eQTL_coloc/seurat/final/coloc_results.xlsx
```

### Current execution logs

```text
processed-data/11_eQTL_coloc/seurat/run_logs/23andMe_2026-07-28/
```

### Frozen old results

```text
processed-data/11_eQTL_coloc/seurat/archive/no23andMe_2026-07-28/
```

### Machine-readable old-versus-new comparisons

```text
processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28/
```

Files there include:

```text
GWAS_input_provenance.tsv
sparse_matched_GWAS_counts.tsv
03_core_invariance.tsv
03_overlap_metrics_long.tsv
03_overlap_metrics_comparison.tsv
03_exact_variant_annotation_changes.tsv
coloc_context_long.tsv
coloc_context_comparison.tsv
coloc_disorder_totals_long.tsv
coloc_gated_strong_changes.tsv
coloc_shared_candidate_PP4_changes.tsv
coloc_category_transitions.tsv
```

## 18. Branch and Filename Policy

Generated tracked workbooks and figures retain canonical names on both
`devel` and `new-gwas`. Git branch identity is the version dimension:

- `devel` has the prior tracked binary content at canonical paths;
- `new-gwas` has regenerated content at the same canonical paths.

Temporary files named `*-new-gwas.<ext>` were created only to preserve local
modified binaries while switching branches. They were moved back to canonical
paths on `new-gwas` before commit and no longer exist.

This is preferable to teaching every generator and downstream consumer a
filename suffix. Side-by-side review uses dated archive/comparison directories.

Important Git limitation: the large raw generated directories are ignored and
are not branch-versioned. Switching branches does not restore ignored raw
outputs. The dated no23 archive is therefore essential for local comparison.

## 19. Observed Downstream Impact

The comparison holds target genotypes, tensorQTL results, DEG definitions,
curated GWAS gene lists, coloc priors, and sensitivity criteria fixed. It
changes the MDD/BD GWAS release and its available QC fields.

Therefore, the deltas should be described as the effect of switching to the
reconstructed integrated GWAS files, not as a perfectly isolated marginal
effect of adding participants.

### 19.1 Genotype-matched GWAS variants

| Trait | `P < 1e-5` no23 | `P < 1e-5` full23 | `P <= 5e-8` no23 | `P <= 5e-8` full23 |
|---|---:|---:|---:|---:|
| BD | 12,660 | 35,186 | 2,894 | 10,470 |
| MDD | 44,555 | 93,770 | 14,053 | 35,668 |

These are target-matched variant rows, not independent loci.

### 19.2 Stable eQTL/DEG rows

| Table | Old rows | New rows | Stable non-GWAS columns equal |
|---|---:|---:|---|
| Significant pairs | 3,435 | 3,435 | yes, tolerance `1e-12` |
| Nominal BH-FDR pairs | 455,143 | 455,143 | yes, tolerance `1e-12` |
| Exploratory GWASx source rows | 3,435 | 3,435 | yes, tolerance `1e-12` |

Curated MDD and BD broad/prioritized gene lists were byte-identical across the
comparison. Thus observed annotation changes came from variant-level GWAS
statistics, not changed eQTL, DEG, or curated gene-list definitions.

### 19.3 Significant eQTL/GWAS overlaps

| Table | Trait | Level | Rows no23 | Rows full23 | Genes no23 | Genes full23 |
|---|---|---|---:|---:|---:|---:|
| Significant pairs | BD | strict | 6 | 11 | 4 | 4 |
| Significant pairs | BD | `P < 1e-5` | 12 | 24 | 7 | 11 |
| Significant pairs | MDD | strict | 8 | 26 | 6 | 12 |
| Significant pairs | MDD | `P < 1e-5` | 24 | 67 | 15 | 39 |
| Nominal BH-FDR | BD | strict | 752 | 2,147 | 13 | 29 |
| Nominal BH-FDR | MDD | strict | 2,593 | 6,466 | 36 | 86 |

### 19.4 Colocalization impact

| Trait | Candidates no23 | Candidates full23 | Gated strong no23 | Gated strong full23 |
|---|---:|---:|---:|---:|
| BD | 584 | 679 | 35 | 22 |
| MDD | 769 | 1,280 | 26 | 46 |

Shared-candidate posterior agreement:

| Trait | Shared candidates | PP4 Pearson | Median absolute PP4 change | Maximum absolute PP4 change | Lead SNP unchanged |
|---|---:|---:|---:|---:|---:|
| BD | 312 | 0.766 | 0.000 | 0.886 | 25.3% |
| MDD | 628 | 0.897 | 0.000 | 0.798 | 47.1% |

Exact gained/lost loci, genes, posterior transitions, and context-level counts
are in the impact report and machine-readable comparison directory. Candidate
rows are screening summaries, not independent GWAS loci.

## 20. Reproduction Guide

### 20.1 Integration test

```bash
cd /dbdata/cdb/gwas-pgs-prs/mbv-prs
bash -n scripts/integrate_23andme_gwas.sh
tests/test_gwas_meta_integration.sh
```

### 20.2 Rebuild integrated GWAS

Use the exact commands in Section 10. Do not use `--resume` after changing an
input, helper script, annotation release, reference, or exclusion list. Use a
new output prefix so old intermediates remain auditable.

### 20.3 Verify installed BCFs

```bash
bcftools index -n /home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.bcf
bcftools index -n /home/gpertea/work/ref/GWAS/MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf

bcftools view -h /home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.bcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP\t%NE\t%SI]\n' \
  /home/gpertea/work/ref/GWAS/BD/bip2024_eur.hg38.bcf | head
```

### 20.4 Downstream harmonization test

```bash
cd /home/gpertea/work/R/spatialDLPFC_mdd_bpd
Rscript -e 'testthat::test_file(
  "code/11_eQTL_coloc/tests/testthat/test-gwas-harmonization.R",
  reporter="summary")'
```

### 20.5 Equivalent downstream rerun commands

The rendered logs prove successful notebook execution, but the exact shell
wrapper/redirection commands were not embedded in every log. Equivalent
commands are:

```bash
cd /home/gpertea/work/R/spatialDLPFC_mdd_bpd

Rscript -e 'rmarkdown::render("code/11_eQTL_coloc/03_eqtl_explore.Rmd")'
Rscript -e 'rmarkdown::render("code/11_eQTL_coloc/03a_nominal_eQTLs.Rmd")'
Rscript -e 'rmarkdown::render("code/11_eQTL_coloc/03b_eQTL_boxplots.Rmd")'
Rscript code/11_eQTL_coloc/03c_GWAS_relaxed_eQTLs.R

Rscript code/11_eQTL_coloc/04_run_coloc.R --disorder MDD --n-cores 12
Rscript code/11_eQTL_coloc/04_run_coloc.R --disorder BD --n-cores 12

Rscript -e 'rmarkdown::render("code/11_eQTL_coloc/05_coloc_explore.Rmd")'
Rscript code/11_eQTL_coloc/07_compare_23andme_impact.R
Rscript -e 'rmarkdown::render("code/11_eQTL_coloc/06_eqtl_coloc_final_tables.Rmd")'
```

Current completion markers cause `04_run_coloc.R` to skip completed contexts.
The runner also refuses partial output sets. A true rerun must first preserve
the current output tree and then use a clean output location/state. Do not
delete current outputs without creating a verified archive.

## 21. Reviewer Checklist

Before accepting these files as a durable analysis release, a new reviewer
should independently verify:

1. Input SHA-256 values against the original delivery and public downloads.
2. BD v7.0 association to v7.2 annotation stability, or obtain v7.0 annotation.
3. Whether an official final DENTIST-filtered BD full-European file can be
   obtained from 23andMe, PGC, or the paper team.
4. Whether an official final MDD full-European file is available for direct
   byte/statistic comparison.
5. Public/23andMe participant independence beyond aggregate cohort counts.
6. Hundreds of allele-oriented sentinel beta, SE, and P values against paper
   supplements, including palindromic variants where possible.
7. Effective-N thresholds and counts directly from unfiltered BCFs.
8. The impact of preserving component-specific INFO/rsq fields and optionally
   applying pre-meta component QC thresholds.
9. The impact of the missing integrated SI rule on downstream variant sets.
10. All liftover-rejected and reference-added records.
11. BCF sample names, headers, record counts, indexes, and checksums.
12. Exact/swapped allele harmonization against the 119-donor PVAR.
13. Current final workbook `gwas_provenance` sheets against BCF hashes.
14. Whether new full `.qs2` coloc objects or explicit SNP-level posterior
    exports are needed for future review.
15. That the coloc dense-cache schema retains `NS`, case fraction is `NC/NS`,
    and a regression test distinguishes it from `NC/NE`.
16. That the installed `coloc` version still ignores `N` and `s` when beta and
    varbeta are supplied, or that outputs are rerun if this behavior changes.
17. That no ignored local work product is mistakenly assumed to be
    branch-versioned or backed up by Git.

## 22. Open Actions

Highest-priority unresolved actions are:

- request the official final full-European BD statistics;
- request the BD DENTIST exclusion list, exact version/parameters, and HRC
  reference if the final file is unavailable;
- request European BD v7.0 annotations or formal v7.0/v7.2 mapping stability;
- request or compare against the official final full-European MDD statistics;
- preserve public INFO and 23andMe avg/min rsq separately in a future
  integration release;
- decide and document component-level imputation QC before interpreting
  differences as participant-only power gains; and
- correct the coloc cache and case-fraction metadata from `NC/NE` to `NC/NS`,
  with a regression test, even though current beta/SE posterior calculations
  are unaffected; and
- rerun coloc with `--write-full-qs2` only if full posterior objects are needed
  and storage cost is justified.

Draft request emails for BD are in:

```text
/dbdata/cdb/gwas-pgs-prs/mbv-prs/BD-DENTIST-issue-v7.0-vs-v7.2.md
```

## 23. Source and Evidence Registry

### Primary literature and methods

- [O'Connell et al. BD paper](https://doi.org/10.1038/s41586-024-08468-9)
- [Adams et al. MDD paper](https://doi.org/10.1016/j.cell.2024.12.002)
- [METAL, Willer et al. 2010](https://doi.org/10.1093/bioinformatics/btq340)
- [GWAS meta-analysis QC, Winkler et al. 2014](https://doi.org/10.1038/nprot.2014.071)
- [GWAS Catalog summary-statistics harmonization](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics)
- [DENTIST](https://github.com/Yves-CHEN/DENTIST)

### Delivered source evidence

```text
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/README_23andMe-update.md
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/7.2-Annotations/v7.2_europe/23andMe_GWAS_Results_v7.2.docx
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/7.2-Annotations/v7.2_europe/23andMe_Platform_Annotations_v7.2.docx
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/Bipolar-Disorder-O_Connell-2025/OConnell_2025_bipolar_european-7.0/bipolar.html
/home/gpertea/work/ref/GWAS/23andMe_MDD_BD/MDD-Adams-2025/Adams_2025_mdd_european-7.2/mdd.html
```

### Integration evidence

```text
/dbdata/cdb/gwas-pgs-prs/mbv-prs/GWAS_23andMe_integration.md
/dbdata/cdb/gwas-pgs-prs/mbv-prs/BD-DENTIST-issue-v7.0-vs-v7.2.md
/home/gpertea/work/ref/GWAS/BD/full_eur_integration/bip2024_eur.integration-report.md
/home/gpertea/work/ref/GWAS/MDD/full_eur_integration/pgc-mdd2025_eur_v3-49-24-11.integration-report.md
/home/gpertea/work/ref/GWAS/BD/full_eur_integration/integration.run.log
/home/gpertea/work/ref/GWAS/MDD/full_eur_integration/integration.run.log
```

### Downstream evidence

```text
code/11_eQTL_coloc/README.md
code/11_eQTL_coloc/Data_Paths.md
code/11_eQTL_coloc/utils.R
code/11_eQTL_coloc/03c_GWAS_relaxed_eQTLs.R
code/11_eQTL_coloc/04_run_coloc.R
code/11_eQTL_coloc/05_coloc_explore.Rmd
code/11_eQTL_coloc/06_eqtl_coloc_final_tables.Rmd
code/11_eQTL_coloc/07_compare_23andme_impact.R
code/11_eQTL_coloc/GWAS-23andMe-MDD-BD-impact-on-tables.md
processed-data/11_eQTL_coloc/seurat/run_logs/23andMe_2026-07-28/
processed-data/11_eQTL_coloc/seurat/archive/no23andMe_2026-07-28/
processed-data/11_eQTL_coloc/seurat/comparison/23andMe_2026-07-28/
```

## 24. Bottom-Line Review Position

The fixed-effect integration method is standard and implemented with robust
allele normalization, explicit sample-size filtering, liftover-aware effect
handling, input hashes, and reusable tests. MDD has good source-release and
sentinel support. BD remains scientifically useful for sensitivity analysis but
must retain its pre-DENTIST and cross-release annotation labels.

The largest review risks are source provenance, missing BD final QC, missing
combined or preserved component imputation-quality fields, and the fact that
downstream differences are not an isolated participant-count contrast. Those
limitations, plus the coloc case-fraction metadata defect documented in
Section 16.4.1, should remain visible in every future report or manuscript use.
