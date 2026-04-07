## Data needs for eQTL and colocalisation analyses
  *   The SpatialExperiment object containing expression data used for Differential Gene Expression (DGE), specifically *before* applying any `filterByExpr` filter.
  *   Genotype data, if available in the project tree.
  *   The Seurat `seurat-pc30` DEG summary files used to compare eQTL findings with DEGs for the same strata.
  *   GWAS summary statistics for MDD and BPD (for eventual colocalization analysis).

## Data locations provided

* **Pivot to Cell-Type (Seurat) Annotations:** The team evaluated Gene Regulatory Network (GRNBoost) outputs and decided that using Seurat label-transfer annotations provided a cleaner biological signal than the initial spatial domains. Consequently, they decided to pivot the eQTL analysis to use these Seurat labels.

Base project path on JHPCE (use ssh `jt` to reach the host with this path)

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
Use only the DEG summaries derived from the same Seurat `seurat-pc30` strata used for the current eQTL analysis:
  - Layer-Adjusted: `processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`
  - Layer-Restricted: `processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv`

These are the relevant DEG files for comparing eQTL findings against DGE in the current strata:
`Astro`, `Inhb`, `L2.3`, `L4`, `L5`, `L6`, `Micro.Vasc`, and `Oligo`.

Older DGE result files derived from other clustering schemes or spatial-domain summaries are not appropriate for per-stratum comparison with the current eQTL results and should be ignored for this purpose.

4. GWAS Summary Statistics (for Colocalization)
Not yet stored/linked in the project directory; check with Shizhong for the proper GWAS data to use.
Details: The team has confirmed they will use European-ancestry-only summary statistics for MDD and Bipolar disorder for the colocalization steps
 They will also evaluate a newer multivariate/cross-disorder psychiatric GWAS from the Psychiatric Genomics Consortium (PGC)
