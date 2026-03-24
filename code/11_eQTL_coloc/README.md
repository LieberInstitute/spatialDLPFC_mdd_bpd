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

Prep workflow:

1. `./stage_required_data.sh` checks/fetches required files from JHPCE (`ssh jt` alias).
2. `./00_get_SNP_PCs.sh` computes genotype PCs (`merged_maf05_pca.eigenvec`).
3. `./01_prep_inputs.R` prepares BED/covariates/exprPC files per `seurat_label`.

Prepared dataset naming convention in `tqtl_in`:

* No `seurat_` prefix
* `all` split has no suffix (e.g., `astro`)
* male/female splits use `_m` and `_f` (e.g., `astro_m`, `astro_f`)
* `micro-vasc` is renamed to `uvasc` (`uvasc`, `uvasc_m`, `uvasc_f`)

Current stratum labels are Seurat-derived (cell-type-like labels):
`Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`/`uvasc`, `Oligo`.

Convenience wrappers:

* `./prepare_tensorqtl_inputs.sh` runs all 3 prep steps above.
* `./prepare_tensorqtl_inputs.sh --check-only` validates preconditions without writing outputs.
* `./02_run_tensorQTL.sh` runs `02a_tensorQTL_cis.py` across all prepared dataset IDs in `tqtl_in` (manifest-driven if available).

Common run patterns:

* `./02_run_tensorQTL.sh --dry-run` to list detected contexts
* `./02_run_tensorQTL.sh --start-from l5_m` to resume from a context
* `./02_run_tensorQTL.sh --only '^(l[2-6]|uvasc)(|_[mf])$'` to restrict contexts
