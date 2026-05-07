This project folder is for code related to eQTL and colocalization analyzes for the MBv project.

## MBv project summary

*   **Study Goal:** The project aims to investigate **sex-specific gene expression changes associated with Major Depressive Disorder (MDD) and Bipolar Disorder (BPD)** across the cortical layers of the dorsolateral prefrontal cortex (dlPFC).
*   **Cohort & Data Collection:** The study utilizes postmortem tissue from **119 adult human brain donors**, evenly split by sex across three diagnostic groups: Neurotypical Controls (NTC, n=40), MDD (n=39), and BPD (n=40).
*   **Technology & Resolution:** Data was generated using the **10x Genomics Visium Spatial RNA-sequencing platform**, capturing spatial transcriptomics (SRT) at a 55um spot resolution.
*   **Data Structure:** The initial unsupervised clustering identified **6 spatial domains** (L1, L2, L3.4, L5, L6, and White Matter). The data from 119 capture areas across 30 slides was pseudobulked by domain and donor, resulting in approximately 690 high-quality samples containing ~14,000 genes for downstream analyses.
*   **Alternative Annotation:** The team explored a Seurat label-transfer approach using single-nucleus RNA-seq to annotate spots by **cell type** (e.g., inhibitory neurons, astrocytes, vasculature) rather than just spatial layers, providing a cleaner biological signal for downstream regulatory networks.

### Differential Gene Expression (DGE) Models

The project DGE analyses use `limma-voom` models on donor-level pseudobulk expression with diagnosis, sex, annotation context, global PC3, age, and `nspots` covariates.

Important terminology note: some information refers to genes with significant `dx*sex interaction`. In the current `code/07_dx_DE` implementation, the final `dx-sex_degs-F-test-t-test.csv` files are generated from sex-stratified diagnosis contrasts such as `F_NTC.MDD`, `M_NTC.MDD`, `F_NTC.BPD`, and `M_NTC.BPD`, followed by an omnibus F-test across the tested contrasts.

A safer wording for these files could be:

- F-test significant dx/sex DGE genes
- genes with sex-stratified diagnosis-associated DGE signal
- broad dx/sex DGE-supported genes

We cannot claim that every gene in these lists has a formal diagnosis-by-sex interaction-only effect unless a contrast matrix explicitly testing only interaction terms is used.

2 related DGE model families are used:

1. **Layer-/context-adjusted (L-A) model:** `~ 0 + condition:sex + context + PC3 + age + nspots`

   This model tests sex-stratified diagnosis contrasts while adjusting for annotation context. Depending on the run, `context` is either a PRECAST smoothed spatial domain or a Seurat label-transfer cluster.

2. **Layer-/context-restricted (L-R) model:** `~ 0 + condition:context:sex + PC3 + age + nspots`

   This model estimates the same diagnosis-by-sex contrast structure within each annotation context, allowing DGE evidence to be localized to a PRECAST spatial domain or a Seurat cluster.

The historical labels `layer-adjusted` and `layer-restricted` are retained in file names. For interpretation in this folder, read `layer` as annotation context:

- PRECAST/smoothed runs use 6 spatial domains: L1, L2, L3.4, L5, L6, WM.
- Seurat runs use 8 cell-type-like clusters: Astro, Inhb, L2.3, L4, L5, L6, Micro.Vasc, Oligo.

### Models to be Considered for eQTL and Colocalization Analyses

The current folder contains code for the eQTL and colocalization stages:

*   **Baseline Model Structure:** Instead of running an "interaction eQTL" model (e.g., `genotype * diagnosis`), which historically yields very few results, the primary strategy will evaluate eQTLs using pseudo-bulking approaches (on sample per donor) to maximize minor allele frequency (MAF > 5%). The model will adjust for diagnosis, donor age, ancestry PCs (assuming mostly European ancestry), and batch effects.
*   **Input Data Pivot:** While initially planned for the 6 spatial domains, the team later decided to run the eQTL analysis on the pseudobulks derived from the **Seurat cell-type labels** (removing the "low UMI" cluster). This aligns the eQTL inputs with the Gene Regulatory Network (GRNBoost/SCENIC) inputs.
*   **Sex-Stratified eQTLs:** Because sex is a major driver of transcriptional differences in this cohort, the eQTL models will be run across three splits: **all donors combined, males only (~60 donors), and females only (~60 donors)**.
*   **Batch Covariates:** The team debated calculating local expression PCs for each group using the `sva` package vs. reusing the global `PC3` variable from the DGE analyses. They concluded that utilizing the existing `PC3` variable as used in DGE acts as a solid shortcut to adjust for technical variation (like slide and sequencing batch) before the eQTL expression PCs capture remaining variance.
*   **Colocalization Analysis:** Once significant spatial/cell-type eQTLs (SNP-gene pairs within a +/- 1 Megabase window) are identified, they will be subjected to colocalization analysis (coloc.abf). This tests if the exact same variant explains both the eQTL expression changes and the genetic risk from GWAS datasets.
*   **GWAS References for Colocalization:** we will compare eQTLs against multiple European-ancestry summary statistics to avoid muddying the signals; independently check **Schizophrenia (SCZ), MDD, and Bipolar GWAS datasets**, as well as a newer **cross-disorder/multivariate psychiatric GWAS**.

## TensorQTL Input Preparation (local `processed-data`)

Canonical local input/output locations for this folder:

* `processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata`
* `processed-data/00_genotypes/plink2/merged_maf05.{pgen,psam,pvar}`
* `processed-data/00_genotypes/plink2/merged_maf05_pca.eigenvec` (generated locally)
* `processed-data/11_eQTL_coloc/tqtl_in` (prepared tensorQTL inputs)
* `processed-data/11_eQTL_coloc/tqtl_out` (tensorQTL outputs)
* `processed-data/11_eQTL_coloc/nspots/tqtl_in` (parallel nspots-adjusted tensorQTL inputs)
* `processed-data/11_eQTL_coloc/nspots/tqtl_out` (parallel nspots-adjusted tensorQTL outputs)

Prep workflow:

1. `./stage_required_data.sh` checks/fetches required files from JHPCE (`ssh jt` alias).
2. `./00_get_SNP_PCs.sh` computes genotype PCs (`merged_maf05_pca.eigenvec`).
3. `./01_prep_inputs.R` prepares BED/covariates/exprPC files per `seurat_label`.

Direct prep entrypoint:

* `Rscript ./01_prep_inputs.R` resolves the built-in defaults, prints them, and prepares outputs.
* `Rscript ./01_prep_inputs.R --help` shows usage plus the exact resolved defaults and exits.
* `Rscript ./01_prep_inputs.R --check-only` performs the full validation/manifest pass without writing outputs.
* `Rscript ./01_prep_inputs.R --dry-run` is an alias for `--check-only`.
* `Rscript ./01_prep_inputs.R --model base` is the default baseline model and writes to `processed-data/11_eQTL_coloc/tqtl_in`.
* `Rscript ./01_prep_inputs.R --model nspots` adds `nspots` to the covariate model, recomputes expression PCs with that expanded design matrix, and writes to `processed-data/11_eQTL_coloc/nspots/tqtl_in` unless `--out-dir` is supplied.
* If staged inputs are missing, `01_prep_inputs.R` reports the missing files and advises running `./stage_required_data.sh` and/or `./00_get_SNP_PCs.sh`.

Covariate models:

* baseline all-donor inputs use `~ DX + sex + age + PC3 + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5`, plus expression PCs.
* baseline male/female inputs use `~ DX + age + PC3 + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5`, plus expression PCs.
* nspots all-donor inputs use `~ DX + sex + age + PC3 + nspots + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5`, plus expression PCs.
* nspots male/female inputs use `~ DX + age + PC3 + nspots + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5`, plus expression PCs.
* `prep_manifest.csv` records the covariate series in `covariate_model`.

Prepared dataset naming convention in `tqtl_in`:

* No `seurat_` prefix
* `all` split has no suffix (e.g., `astro`)
* male/female splits use `_m` and `_f` (e.g., `astro_m`, `astro_f`)
* `micro-vasc` is renamed to `uvasc` (`uvasc`, `uvasc_m`, `uvasc_f`)
* The nspots series keeps the same dataset IDs; the directory path distinguishes baseline from nspots-adjusted inputs and outputs.

Current stratum labels are Seurat-derived (cell-type-like labels):
`Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`/`uvasc`, `Oligo`.

Convenience wrappers:

* `./prepare_tensorqtl_inputs.sh` is deprecated but still available as a shim that runs all 3 prep steps above.
* `./prepare_tensorqtl_inputs.sh --check-only` is deprecated but still validates preconditions without writing outputs.
* `./02_run_tensorQTL.sh` runs `02a_tensorQTL_cis.py` across all prepared dataset IDs in `tqtl_in` (manifest-driven if available).
* `./02_run_tensorQTL.sh --analysis nspots` runs the same mapper against `processed-data/11_eQTL_coloc/nspots/tqtl_in` and writes to `processed-data/11_eQTL_coloc/nspots/tqtl_out`.
* `./02_run_tensorQTL.sh --input-dir PATH --output-dir PATH` can run any explicitly supplied prepared input/output pair.

Common run patterns:

* `./02_run_tensorQTL.sh --dry-run` to list detected contexts
* `./02_run_tensorQTL.sh --start-from l5_m` to resume from a context
* `./02_run_tensorQTL.sh --only '^(l[2-6]|uvasc)(|_[mf])$'` to restrict contexts
* `Rscript ./01_prep_inputs.R --model nspots --check-only` to validate the nspots design without writing files
* `Rscript ./01_prep_inputs.R --model nspots` to prepare all nspots-adjusted inputs
* `./02_run_tensorQTL.sh --analysis nspots --dry-run` to list nspots-adjusted tensorQTL commands
* `./02_run_tensorQTL.sh --analysis nspots --only '^astro$'` to run one nspots-adjusted dataset

### Custom-cluster tensorQTL inputs

The custom-cluster workflow is an additional analysis series and does not reuse or overwrite the existing Seurat `tqtl_in` or nspots `nspots/tqtl_in` inputs.

Custom prep uses the JHPCE-derived object:

* `processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata`

Primary custom model:

* grouping column: `custom_cluster`
* split: `all`
* covariates: `DX + sex + age + PC3 + nspots + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5`, plus expression PCs
* output input directory: `processed-data/11_eQTL_coloc/custom_cluster_nspots/tqtl_in`
* tensorQTL output directory: `processed-data/11_eQTL_coloc/custom_cluster_nspots/tqtl_out`

Prepare custom inputs:

```bash
./stage_required_data.sh
Rscript ./01_prep_inputs.R \
  --spe-file ../../processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata \
  --cluster-col custom_cluster \
  --analysis-label custom_cluster_nspots \
  --model nspots \
  --splits all \
  --out-dir ../../processed-data/11_eQTL_coloc/custom_cluster_nspots/tqtl_in
```

Run custom tensorQTL:

```bash
./02_run_tensorQTL.sh --analysis custom_cluster_nspots
```

Dry-run only:

```bash
./02_run_tensorQTL.sh --analysis custom_cluster_nspots --dry-run
```

## DGE Comparison Inputs

Use the final DEG summary CSVs as the source of DGE significance, with DEG status defined by the F-test BH-adjusted p-value (`adj.P.Val < 0.05`) from the summary files.

Do not add a global `n_ttest_sig > 0` filter when reconstructing Jacqui's broad DEG universe. The post-hoc t-test columns are useful for localizing and interpreting effects, but Jacqui's recommended broad list is based on F-test significance.

Recommended language:

- Use: "F-test significant dx/sex DGE genes"
- Use: "sex-stratified diagnosis-associated DGE genes"
- Use: "broad PRECAST-or-Seurat DGE-supported genes"
- Avoid: "all genes with a formal dx-by-sex interaction effect"

The shorthand `dx*sex DEG` is useful for matching existing project file names and Slack language, but it should not be overinterpreted as proving a formal interaction-only effect for every gene.

Global DEG support for this project comes from the union of all four author-recommended DEG summaries in `processed-data/07_dx_DE`:

- PRECAST L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- PRECAST L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
- Seurat L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- Seurat L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

The author-provided union list at `raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt` is useful as a validation artifact for the reconstructed four-file union by `gene_name`.

The tensorQTL results in this folder are Seurat-context eQTLs. Therefore, DEG/eQTL overlaps have two distinct interpretations:

1. **Global DEG support**

   Ask whether an eQTL eGene is present in the broad PRECAST-or-Seurat DEG universe. This is the closest match to Jacqui's Slack guidance to "paint with a broad brush."

2. **Seurat context support**

   Ask whether an eQTL eGene is DEG-supported in the same Seurat context used for tensorQTL. This is more specific, but it uses only the Seurat DEG files because PRECAST domains do not map one-to-one onto the 8 Seurat eQTL strata.

Use `gene_id` as the primary overlap key against tensorQTL phenotypes. Keep `gene_name` for reporting and validation against the author text list.

### Additional DGE list notes/info

For now the goal is to paint with a broad brush, and so we would like to use all of the genes with sig. dx*sex interaction (F-test adj, p<.05) in either the PRECAST (smoothed) annotation (6 domains) or the Seurat annotations (8 clusters).

  - The layer-adjusted model results are highly concordant (D-E in attached image)
  - The layer-restricted model results are mostly concordant (F-G in attached image)
  - The pseudobulk samples you are working with share characteristics of both annotation sets. Like the PRECAST domains there are no spots that were annotated to the low UMI cluster. Like the Seurat model results, the spots are aggregated based on cell-type labels in to 8 clusters. For the Seurat DE models, the low UMI cluster spots were included. So the pseudobulk data you have is identical to neither.
  - We used all of the genes with sig. dx*sex interaction in either model when we constructed our DEG modules.

A txt list of all of these genes is present here:
 `raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt`
That list is a compilation (union) of the two files you specified plus these two files:
  - Layer-adjusted domains: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - Layer-restricted domains: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`

### DEG views used for eQTL overlap summaries

The overlap code can report 4 DEG views.
For preliminary reporting, prioritize the first two views:

1. `broad_interaction`
2. `context_localized`

The sex-aware views are useful secondary or backup analyses for male/female tensorQTL splits:

3. `sex_specific`
4. `context_and_sex_specific`

File mapping:

| DEG view | DEG files used | Extra DEG filter after `adj.P.Val < 0.05` | eQTL matching |
| --- | --- | --- | --- |
| `broad_interaction` | PRECAST L-A, PRECAST L-R, Seurat L-A, Seurat L-R | none | gene only |
| `context_localized` | Seurat L-R only | `n_ttest_sig_<Seurat context> > 0` | gene and Seurat context |
| `sex_specific` | PRECAST L-A, PRECAST L-R, Seurat L-A, Seurat L-R | any same-sex post-hoc t-test flag: `F_*_ttest` for female or `M_*_ttest` for male | gene and sex |
| `context_and_sex_specific` | Seurat L-R only | same-context, same-sex post-hoc t-test flag, for example `Astro_F_NTC.MDD_ttest` | gene, Seurat context, and sex |

All four views start from F-test significant rows in the named files (`adj.P.Val < 0.05`). The post-hoc t-test columns only localize the already F-test significant DEG support by Seurat context and/or sex.

#### 1. `broad_interaction`: inclusive project-level DGE support

Definition: Union of F-test significant genes from all 4 recommended DEG files:
  - PRECAST L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - PRECAST L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - Seurat L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
  - Seurat L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

Matching to eQTLs:
- Not context-matched.
- Not sex-matched.
- Gene-level overlap only.

This is a primary broad screening view, "painting with a broad brush".
- Answers: is this eQTL eGene in the project-level dx/sex DGE universe?

#### 2. `context_localized`: same Seurat context support

Definition:
- Source file: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- Filter: Seurat L-R F-test significant genes localized to the same Seurat context via `n_ttest_sig_<context> > 0`.

Matching to eQTLs: context-matched, not sex-matched.

This is a primary context-aware view for all-donor Seurat eQTLs, showing eGenes overlapping Seurat L-R DEGs localized to the same cluster.
- Answers: is this eQTL eGene DEG-supported in the same Seurat cluster?

Note: this view excludes PRECAST-only DGE support because PRECAST domains do not map cleanly onto the Seurat eQTL contexts.

#### 3. `sex_specific`: same-sex support without context matching

Definition: Genes from all four DEG tables with a post-hoc sex-specific t-test flag matching the eQTL split:
  - PRECAST L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - PRECAST L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - Seurat L-A: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
  - Seurat L-R: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
  - female eQTLs use `F_*_ttest`
  - male eQTLs use `M_*_ttest`

Matching to eQTLs: Sex-matched, not context-matched. Secondary view for sex-stratified eQTLs.
- Answers: does a male or female eQTL eGene overlap any same-sex DGE-supported gene, regardless of annotation context?

This is not a context-localized result, so it is broader than the Seurat eQTL context.

#### 4. `context_and_sex_specific`: strict same Seurat context and same sex support

Definition:
- Source file: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
- Filter: Seurat L-R genes with both matching context and matching sex-specific t-test evidence.
- Example: for female Astro eQTLs, use columns such as `Astro_F_NTC.MDD_ttest`.

Matching to eQTLs: Context-matched, sex-matched.
- Strictest interpretive view for sex-stratified Seurat eQTLs.
- Answers: does a sex-stratified eQTL eGene overlap a DEG in the same Seurat cluster and same sex?

This view is expected to have smaller counts because it requires both same-context and same-sex DEG support.


## R coding agent instructions
  - use single line comments starting with '## ' and lower case, to briefly comment/explain non-trivial code blocks generated
  - use here::i_am('.git/HEAD') in R/Rmd to anchor the project base folder, make all project paths relative to it
  - prefer base R and data.table over dplyr for data manipulation, reshaping, filtering etc.; avoid local/global variable name conflicts/clash with data.table columns which can lead to serious silent logic bugs with data.table notation
