This project folder is for code related to eQTL and colocalization analyzes for the MBv project.

## MBv project summary

*   **Study Goal:** The project aims to investigate **sex-specific gene expression changes associated with Major Depressive Disorder (MDD) and Bipolar Disorder (BPD)** across the cortical layers of the dorsolateral prefrontal cortex (dlPFC).
*   **Cohort & Data Collection:** The study utilizes postmortem tissue from **119 adult human brain donors**, evenly split by sex across three diagnostic groups: Neurotypical Controls (NTC, n=40), MDD (n=39), and BPD (n=40).
*   **Technology & Resolution:** Data was generated using the **10x Genomics Visium Spatial RNA-sequencing platform**, capturing spatial transcriptomics (SRT) at a 55um spot resolution.
*   **Data Structure:** The initial unsupervised clustering identified **6 spatial domains** (L1, L2, L3.4, L5, L6, and White Matter). The data from 119 capture areas across 30 slides was pseudobulked by domain and donor, resulting in approximately 690 high-quality samples containing ~14,000 genes for downstream analyses.
*   **Alternative Annotation:** The team explored a Seurat label-transfer approach using single-nucleus RNA-seq to annotate spots by **cell type** (e.g., inhibitory neurons, astrocytes, vasculature) rather than just spatial layers, providing a cleaner biological signal for downstream regulatory networks.

### Differential Gene Expression (DGE) Models

The team utilized a `limma` framework to leverage the interaction between diagnosis and sex (`dx*sex`), maximizing the statistical power of the 119 donors. Two primary models were designed, both utilizing covariates for technical batch effects (global PC3), age, and the number of spots per pseudobulk sample (`nspots`):
1.  **Layer-Adjusted (L-A) Model:** `~ dx*sex + domain + PC3 + age + nspots`. This model tests the overall influence of diagnosis and sex on gene expression while statistically adjusting for the differences between spatial layers.
2.  **Layer-Restricted (L-R) Model:** `~ dx*sex*domain + PC3 + age + nspots`. This model directly probes the interaction between diagnosis, sex, and spatial domain, allowing the team to identify gene expression changes restricted to specific layers (e.g., White Matter alterations specifically in females).

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

## DGE Comparison Inputs

Use the final DEG summary CSVs as the source of significance, with DEG status defined by the F-test BH-adjusted p-value (`adj.P.Val < 0.05`) from the summary files. Do not add a global `n_ttest_sig > 0` filter.

Global DEG support for this project comes from the union of all four author-recommended DEG summaries in `processed-data/07_dx_DE`:

* `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
* `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
* `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
* `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

The author-provided union list at `raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt` is useful as a validation artifact for the reconstructed four-file union by `gene_name`.

Per-context eQTL overlap remains anchored to the current Seurat strata:
`Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`, and `Oligo`.

Per-context DEG sets should be defined as:

* all Seurat `layer-adjusted` F-test-significant genes, plus
* Seurat `layer-restricted` F-test-significant genes localized to that context via `n_ttest_sig_<context> > 0`

This keeps PRECAST/domain DEG support global while avoiding a forced one-to-one mapping from PRECAST domains onto the 8 Seurat contexts used by the current eQTL workflow.

Use `gene_id` as the primary overlap key against tensorQTL phenotypes. Keep `gene_name` for reporting and validation against the author text list.

The default `n_context_DEGs` overlap summaries use this broad per-context reference, so all/male/female eQTL splits share the same DEG reference set and several contexts are close to the 523 Seurat layer-adjusted genes. Additional DEG-view summaries in `03_eqtl_explore.Rmd` separate stricter interpretations: `broad_interaction`, `context_localized`, `sex_specific`, and `context_and_sex_specific`. Context-aware views use only the eight Seurat eQTL contexts above, not PRECAST/smoothed domains.

DEG-view matching to eQTL strata:

* `broad_interaction` is the main Jacqui-recommended broad support set. It is not context- or sex-matched to a specific eQTL stratum; it is best for inclusive screening.
* `context_localized` is the best context-matched view for all-donor eQTLs. Its DEG provenance is only the Seurat layer-restricted table, filtered to `adj.P.Val < 0.05` and `n_ttest_sig_<context> > 0`. It is not derived from male/female DEG lists, so the same context set is used for all, male, and female eQTL splits.
* `sex_specific` is sex-matched but not context-matched. It uses all four DEG tables and any post-hoc sex-specific t-test flag matching `F_*_ttest == "padj<.05"` or `M_*_ttest == "padj<.05"`. In all-donor eQTL summaries, female and male DEG sets are reported as separate `deg_sex` rows; those are the same female and male DEG sets used for the female and male eQTL splits.
* `context_and_sex_specific` is the closest match for sex-stratified eQTLs because it matches both the Seurat eQTL context and the eQTL sex split. Its DEG provenance is only the Seurat layer-restricted table and context-sex columns such as `Astro_F_NTC.MDD_ttest == "padj<.05"`. In all-donor eQTL summaries, female and male context-sex DEG sets are reported separately rather than merged.

### Additional DGE list note

[Jacqui on Slack:]

I think right now the goal is to paint with a broad brush, and so we would like to use all of the genes with sig. dx*sex interaction (F-test adj, p<.05) in either the PRECAST (smoothed) annotation (6 domains) or the Seurat annotations (8 clusters). Here are a few reasons:

  - The layer-adjusted model results are highly concordant (D-E in attached image)
  - The layer-restricted model results are mostly concordant (F-G in attached image)
  - The pseudobulk samples you are working with share characteristics of both annotation sets. Like the PRECAST domains there are no spots that were annotated to the low UMI cluster. Like the Seurat model results, the spots are aggregated based on cell-type labels in to 8 clusters. For the Seurat DE models, the low UMI cluster spots were included. So the pseudobulk data you have is identical to neither.
  - We used all of the genes with sig. dx*sex interaction in either model when we constructed our DEG modules.


A txt list of all of these genes is present here:
 `raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt`
That list is a compilation (union) of the two files you specified plus these two files:
  - Layer-adjusted domains: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - Layer-restricted domains: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`

## R coding agent instructions
  - use single line comments starting with '## ' and lower case, to briefly comment/explain non-trivial code blocks generated
  - use here::i_am('.git/HEAD') in R/Rmd to anchor the project base folder, make all project paths relative to it
  - prefer base R and data.table over dplyr for data manipulation, reshaping, filtering etc.; avoid local/global variable name conflicts/clash with data.table columns which can lead to serious silent logic bugs with data.table notation
