## Data needs for eQTL and colocalisation analyses
  *   The SpatialExperiment object containing expression data used for Differential Gene Expression (DGE), specifically *before* applying any `filterByExpr` filter.
  *   Genotype data, if available in the project tree.
  *   The final DEG summary files and author DEG union list used to compare eQTL findings with DEGs.
  *   GWAS summary statistics for MDD and BPD (for eventual colocalization analysis).

## Data locations provided

* **Pivot to Cell-Type (Seurat) Annotations:** The team evaluated Gene Regulatory Network (GRNBoost) outputs and decided that using Seurat label-transfer annotations provided a cleaner biological signal than the initial spatial domains. Consequently, they decided to pivot the eQTL analysis to use these Seurat labels.

Base project path on JHPCE (use ssh `jh` to reach the host with this path)

`/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/`

The team decided to pivot the downstream analyses—including the eQTL and Gene Regulatory Network (GRN) models—to use the Seurat label transfer annotations rather than the original spatial domains.
Here is the current summary of the exact data locations and variables required for Geo's tensorQTL and colocalization runs:

1. Revised Pseudobulk Expression Object (primary input expression data):
 `processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata`
Details: This updated SpatialExperiment object uses the Seurat cell-type label annotations and completely removes the "low UMI" cluster spots so that the eQTL inputs perfectly match the GRN/SCENIC analysis
Key Variables: The cluster information used to stratify the eQTL runs is stored in colData under the variable seurat_label
The model design remains identical to previous plans, so it still relies on the PC3 variable to adjust for batch effects

2. Genotype Data
These files are symlinked/copied in `processed-data/00_genotypes` folder;
MAF>0.05 files should be used (plink2 prefix `plink2/merged_maf05`)

3. Differential Gene Expression (DEG) Lists for eQTL overlap assessment
DEG significance should be taken directly from the final summary CSVs using `adj.P.Val < 0.05` from the F-test summaries. Do not add a global `n_ttest_sig > 0` filter.

Global DEG support should be reconstructed from the union of these four DEG summaries:
  - `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv`
  - `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
  - `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

The author-provided DEG union list is:
  - `raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt`

Use that text file as a validation target for the reconstructed four-file union by `gene_name`.

Per-context overlap with the current eQTL strata should remain Seurat-context-specific:
`Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`, and `Oligo`.

Per-context DEG sets should be defined as:
  - all Seurat layer-adjusted F-test-significant genes, plus
  - Seurat layer-restricted F-test-significant genes with `n_ttest_sig_<context> > 0`

This keeps PRECAST/domain DEG support global, while localization to current eQTL contexts stays tied to the Seurat labels actually used in tensorQTL.

Use `gene_id` as the primary overlap key against tensorQTL phenotype IDs, with `gene_name` retained for reporting and validation.

4. GWAS Summary Statistics (for Colocalization)
Not yet stored/linked in the project directory; check with Shizhong for the proper GWAS data to use.
Details: The team has confirmed they will use European-ancestry-only summary statistics for MDD and Bipolar disorder for the colocalization steps
 They will also evaluate a newer multivariate/cross-disorder psychiatric GWAS from the Psychiatric Genomics Consortium (PGC)
