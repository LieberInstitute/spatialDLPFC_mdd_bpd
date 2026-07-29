# GWAS Summary Statistics on GRCh38

This directory contains psychiatric GWAS summary statistics in GWAS-VCF BCF
format on GRCh38 and matched PGS loading outputs from the BCFtools `+pgs`
plugin.

The source summary statistics were normalized with `bcftools +munge`, lifted from GRCh37/hg19 to GRCh38 with `bcftools +liftover`, sorted, and indexed. PGS loadings were computed from the converted GWAS BCFs with the EUR LDGM file and `bcftools +pgs`.

All main outputs are BCF files with matching `.bcf.csi` indexes. They use GRCh38 contig names with `chr` prefixes.

## Current Analysis Files

| Disorder | Converted GWAS BCF | Sample | Records | PGS BCF | PGS sample | PGS records |
|---|---|---:|---:|---|---:|---:|
| BD | `BD/bip2024_eur.hg38.bcf` | `BD_2024_FULL_EUR` | 6394788 | `BD/bip2024_eur.hg38.pgs.b5e-8.bcf` | `BD_2024_FULL_EUR_pgs_a0.5_b5e-08` | 5227123 |
| MDD | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf` | `MDD_2025_FULL_EUR` | 6656222 | `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.pgs.b2e-8.bcf` | `MDD_2025_FULL_EUR_pgs_a0.5_b2e-08` | 5339102 |
| SCZD | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf` | `SCZ_2022.EUR` | 7658487 | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf` | `SCZ_2022.EUR_pgs_a0.5_b2e-07` | 6076466 |

The BD and MDD files integrate the public no-23andMe European summary
statistics with the supplied European 23andMe component by fixed-effect
inverse-variance meta-analysis. The BD reconstruction is pre-DENTIST and must
not be called the exact final paper GWAS. The public no-23andMe files remain in
the same directories for sensitivity comparisons.

The current integrated BCF checksums are:

```text
1d502351659d81aa503da64779b7cc9e80c6819bca9e2ab5d5ac6efc3bb531be  BD/bip2024_eur.hg38.bcf
a8b30df37b032097920ded06697aa60851f0078e681f1a595b6d5f3ec70347f3  MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf
```

## Required Files for This Project

Current eQTL/GWAS overlap and `loadGWAS()` cache rebuilding require:

- `BD/bip2024_eur.hg38.bcf`
- `BD/bip2024_eur.hg38.bcf.csi`
- `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf`
- `MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf.csi`
- `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf`
- `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.bcf.csi`
- standardized and source GWAS gene-list TSVs in each disorder folder.

`coloc.abf` uses the same full `.hg38.bcf` files and indexes because they
contain genome-wide `ES`, `SE`, `LP`, `NE`, `NS`, and `NC` fields on GRCh38.
The reconstructed integrated BD and MDD BCFs have no combined `SI`; this is
represented as missing, not imputed or borrowed from either component. The
project therefore applies `SI >= 0.8` to SCZD but no post-integration SI filter
to BD or MDD.

The PGS outputs are not required for current eQTL/GWAS overlap or future `coloc.abf`:

- `*.hg38.pgs.*.bcf`
- `*.hg38.pgs.*.bcf.csi`
- `*.hg38.pgs.*.log`

Original downloaded summary-statistic archives and `colheaders.tsv` are provenance/rebuild inputs and are optional if the full `.hg38.bcf` files are already staged.

## Source GWAS Data

SCZD uses the public European PGC summary statistics. BD and MDD use the
integrated products described above; their public no-23andMe components remain
available from the original paper repositories.

| Disorder | Local source file used | Upstream file | Paper and DOI | Download location |
|---|---|---|---|---|
| BD | `BD/bip2024_eur_no23andMe.gz` plus approved 23andMe European component | `bip2024_eur_no23andMe.gz` and O'Connell 2025 23andMe delivery | O'Connell KS, Koromina M, van der Veen T, et al. Genomics yields biological and phenotypic insights into bipolar disorder. Nature. doi:10.1038/s41586-024-08468-9 | Figshare record https://doi.org/10.6084/m9.figshare.27216117 |
| MDD | `MDD/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz` plus approved 23andMe European component | public no-23andMe file and Adams 2025 23andMe delivery | Adams MJ, Streit F, Meng X, Awasthi S, et al. Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies. Cell. doi:10.1016/j.cell.2024.12.002 | Figshare record https://doi.org/10.6084/m9.figshare.27061255 |
| SCZD | `SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz` | `PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz` | Trubetskoy V, Pardinas AF, Qi T, et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. Nature. doi:10.1038/s41586-022-04434-5 | Figshare record https://doi.org/10.6084/m9.figshare.19426775; direct file https://ndownloader.figshare.com/files/34517828 |

## VCF Columns and IDs

Each BCF uses standard VCF columns:

| Column | Meaning |
|---|---|
| `CHROM` | GRCh38 chromosome or contig, usually `chr1` to `chr22` for these autosomal outputs |
| `POS` | 1-based GRCh38 position |
| `ID` | source variant identifier from the GWAS input |
| `REF` | GRCh38 reference allele after liftover |
| `ALT` | alternate/effect allele used by GWAS-VCF FORMAT fields |
| `QUAL` | unset in these files |
| `FILTER` | site filter status |
| `INFO` | liftover metadata, when applicable |
| `FORMAT` | per-sample GWAS or PGS fields |
| sample column | one synthetic sample holding summary statistics |

If a dbSNP rsID is present, it is in VCF column 3, `ID`. Not every source identifier is an rsID: MDD and SCZD include some `chr:pos_ref_alt` style IDs. Current ID coverage:

| File type | BD rsIDs | MDD rsIDs | SCZD rsIDs |
|---|---:|---:|---:|
| Converted GWAS BCF | 6360316 / 6394788 | 6645867 / 6656222 | 7637489 / 7658487 |
| PGS BCF | 5223603 / 5227123 | 5331173 / 5339102 | 6075044 / 6076466 |

For harmonization with genotype VCF/BCF files, match on GRCh38 `CHROM`, `POS`, `REF`, and `ALT`. Treat `ES` as relative to `ALT`.

## Converted GWAS BCF FORMAT Fields

The converted `.hg38.bcf` files store GWAS summary statistics in FORMAT fields. Field values are per alternate allele (`Number=A`).

| FORMAT field | Meaning |
|---|---|
| `NS` | variant-specific number of samples or individuals with called genotypes |
| `SI` | imputation accuracy score; unavailable after the BD/MDD integration |
| `NC` | variant-specific number of cases |
| `ES` | effect size estimate relative to `ALT`; BETA sources stay beta, OR sources are converted by `+munge` to log effect |
| `SE` | standard error of `ES` |
| `LP` | `-log10(P)` for the effect estimate |
| `NE` | variant-specific effective sample size |
| `I2` | Cochran heterogeneity I squared |
| `CQ` | Cochran Q `-log10(P)` |
| `ED` | effect direction across studies |

Per-disorder FORMAT layouts:

| Disorder | FORMAT fields in converted GWAS BCF |
|---|---|
| BD | `NS:NC:ES:SE:LP:NE:I2:CQ:ED`; `SI` queries return missing |
| MDD | `NS:NC:ES:SE:LP:NE:I2:CQ:ED`; `SI` queries return missing |
| SCZD | `NS:SI:NC:ES:SE:LP:NE` |

## PGS BCF FORMAT Fields

The `.hg38.pgs.*.bcf` files are the `bcftools +pgs` outputs. They contain:

| FORMAT field | Meaning |
|---|---|
| `ES` | GraphPred PGS loading/effect score relative to `ALT` |

The PGS outputs were created with EUR LDGM and these options:

| Disorder | `--beta-cov` | `--max-alpha-hat2` | Input exclusion |
|---|---:|---:|---|
| BD | `5e-8` | `0.001` | `FILTER="IFFY"` |
| MDD | `2e-8` | `0.0005` | `FILTER="IFFY"` |
| SCZD | `2e-7` | `0.002` | `FILTER="IFFY"` |

## INFO and FILTER Fields

The converted GWAS and PGS BCFs share liftover INFO definitions:

| INFO field | Meaning |
|---|---|
| `FLIP` | allele was strand-flipped during liftover |
| `SWAP` | alternate allele became reference during liftover; `-1` means a new reference allele was added |

Most records have `INFO=.`. Records with liftover changes may have `FLIP`, `SWAP=1`, `SWAP=-1`, or both.

The shared FILTER definitions are:

| FILTER | Meaning |
|---|---|
| `IFFY` | reference allele could not be determined |
| `REF_MISMATCH` | reference does not match any allele |

The inspected outputs currently have all records unfiltered as `FILTER=.`.

## colheaders.tsv

`colheaders.tsv` is the shared two-column header map used by `bcftools +munge` to translate source summary-statistic headers into canonical GWAS-VCF inputs. It is not disorder-specific.

Key mappings:

| Source headers | Canonical field |
|---|---|
| `SNP`, `ID` | `SNP` |
| `CHR`, `CHROM`, `#CHROM` | `CHR` |
| `BP`, `POS` | `BP` |
| `A1`, `EA` | `A1` |
| `A2`, `NEA` | `A2` |
| `OR` | `OR` |
| `BETA` | `BETA` |
| `SE` | `SE` |
| `P`, `PVAL` | `P` |
| `INFO`, `IMPINFO` | `INFO` |
| `NCAS`, `Nca` | `N_CAS` |
| `NCON`, `Nco` | `N_CON` |
| `NEFF` | `NEFF` |
| `NEFFDIV2`, `Neff_half` | `NEFFDIV2` |
| `HETI`, `HetISqt` | `HET_I2` |
| `HETPVAL`, `HetPVa` | `HET_P` |
| `Direction` | `DIRE` |
| `HRC_FRQ_A1` | `FRQ` |

`FCAS` and `FCON` are intentionally not mapped because they are case and control allele frequencies, not one overall effect allele frequency.

## Inspection Commands

Useful checks:

```bash
bcftools view -h BD/bip2024_eur.hg38.bcf
bcftools query -l MDD/pgc-mdd2025_eur_v3-49-24-11.hg38.bcf
bcftools index -n SCZD/PGC3_SCZ_wave3.european.autosome.public.v3.hg38.pgs.b2e-7.bcf
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES]\n' BD/bip2024_eur.hg38.pgs.b5e-8.bcf | head
```
