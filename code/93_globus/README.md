### Analysis Part 1: Spot-level processing
| File name | Description |
| :----- | :----- |
| MBv_n119_features.csv.gz | 36,601 genes from `GRCh38 2020-A`, annotated with number of spots with non-zero expression. For convenience, includes columns for important gene filters determined by downstream analyses. |
| MBv_n119_observations.csv.gz |  594,044 spots from 119 donors containing: <ul><li>Donor metadata: sample_id (consisting of visium slide and position), brnum, MBv_sample, diagnosis, biological sex, age</li><li>Sequencing batch: seq</li><li>Spot information: unique barcode key, array_row, array_col, in_tissue</li><li>QC metrics: sum_umi (library size), sum_gene (detected genes), expr_chrM_ratio (mito. fraction), remove_spots (QC cutoff)</li><li>Domain annotation: domain_sp, domain_ct</li></ul> |
| MBv_n119_spatial-coords.csv.gz | Space Ranger spatial output for 594,044 spots from 119 donors. |
| MBv_n119_raw-counts.csv.gz | Raw counts for 36,601 genes (rows) and 594,044 spots (columns). |

Files with prefix `MBv_n119_` can be used to construct the full unfiltered dataset and to replicate our analysis following the code on our [Github repo](https://research.libd.org/spatialDLPFC_mdd_bpd/). These files are also available as supplementary data associated with accession number [GSE333974](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333974).

As an alternative to re-running all of the code, `MBv_n119_observations.csv` includes spot-level metadata for key steps of the processing pipeline. 
Users can filter to the 535,248 spots that passed QC (indicated by 'remove_spots') and utilize the domain<sub>SP</sub> and domain<sub>CT</sub> annotation labels provided for these spots to explore spot-level data or perform custom pseudobulking.

### Analysis Part 2: Pseudobulk for differential expression (DE)
| File name | Description |
| :----- | :----- |
| MBv_pseudobulk_domain-sp.rds | Fully processed (QC filtered and `logcounts` computed) pseudobulk `SummarizedExperiment` object for domain<sub>SP</sub> annotations (21,080 genes, 690 samples). Includes `colData` with pseudobulk QC metrics, `rowData` with columns labeling the 1663 SVGs and the subset of 13,844 genes input for DE analysis. Includes PCA in `reducedDim` slot. |
 MBv_pseudobulk_domain-ct.rds | Fully processed (QC filtered and `logcounts` computed) pseudobulk `SummarizedExperiment` object for domain<sub>CT</sub> annotations (21,175 genes, 942 samples). Includes `colData` with pseudobulk QC metrics, `rowData` with columns labeling the 1663 SVGs and the subset of 13,119 genes input for DE analysis. Includes PCA in `reducedDim` slot. |

Users can directly access the pseudobulk SummarizedExperiment objects for domain<sub>SP</sub> or domain<sub>CT</sub> annotations with the files `MBv_pseudobulk_domain-sp.rds` or `MBv_pseudobulk_domain-ct.rds`, respectively.

The objects reflect the pseudobulk-level QC filters and include the `logcount` values, as well as the PCA results.
The pseudobulk objects can be used to replicate the DE models implemented in this study, and the gene filters used for DE can be found in the `rowData` slot or in `MBv_n119_features.csv`. 

### Analysis Part 3: GRNBOOST2, Pseudobulk for eQTL 
| File name | Description |
| :----- | :----- |
| MBv_pseudobulk_domain-ct-eQTL.rds | Fully processed (QC filtered and `logcounts` computed) pseudobulk `SummarizedExperiment` object for domain<sub>CT</sub> annotations, after removing spots still annotated to domain<sub>SP</sub> low UMI cluster after smoothing (21,174 genes, 943 samples). Includes `colData` with pseudobulk QC metrics, `rowData` with columns labeling the 1663 SVGs and the subset of 13,162 genes input for `GRNBOOST2` analysis. Includes PCA in `reducedDim` slot. |
| MBv_GRNBOOST2_F-test-padj-05_edge-list.csv | `GRNBOOST2` output used to generate DEG modules from 119 donors. Predictors used: 817 genes with F-test adj. p<.05 in any DE model. Includes columns of predictor, target, importance, regulation direction, and rho. |
| MBv_GRNBOOST2_human-TF_edge-list.csv | `GRNBOOST2` output used to generate transcription factor (TF) regulons from 109 donors. Predictors used: 1,892 human TFs (`SCENIC` default). Includes columns of TF predictor, target, importance, regulation direction, and rho. |
| MBv_n119_DEG-modules_AUCell.csv.gz | `AUCell` output providing spot-level estimation of DEG module enrichment (513,200 spots, 15 modules). |
| MBv_n109_TF-regulons_AUCell.csv.gz | `AUCell` output providing spot-level estimation of TF regulon enrichment (469,396 spots, 18 regulons). |

Module, regulon, and eQTL analysis all utilized a combination of domain<sub>SP</sub> and domain<sub>CT</sub> information. 
Prior to these analyses and after QC, spots annotated to the domain<sub>SP</sub> `low UMI` cluster (after smoothing) were removed, reducing the number of spots to 513,200.
For eQTL analysis, these spots were then pseudobulked according to domain<sub>CT</sub> annotation, and underwent pseudobulk-level QC, calculation of `logcounts`, and gene filtering. 
The resulting pseudobulk `SummarizedExperiment` object (21,174 genes, 943 samples) can be accessed with the `MBv_pseudobulk_domain-ct-eQTL.rds` file.

The raw counts matrix used as input for `GRNBOOST2` (module and regulon analysis), can be generated from the spot-level raw data files ('Part 1'). 
Simply use obsercation columns 'remove_spots' to isolate the spots that passed QC and use 'domain_sp' to remove to spots assigned to the low UMI cluster ('dropped.lowUMI') to identify the 513,200 spots used as input for module analysis. For regulon analysis, use the additional spot-level filter provided in 'SCENIC_n109'. 

The 13,162 genes input to `GRNBOOST2` resulted from the `edgeR` filters from the `MBv_pseudobulk_domain-ct-eQTL.rds` file and can be found in the `rowData` slot of this object or in the `MBv_n119_features.csv` column 'grnboost2'. 

Due to the stochastic nature of `GRNBOOST2` output resulting from the use of random number seeds, every time this analysis is run (with identical input and parameters) it will produce slightly different outcomes. 
Therefore we have provided for users the unfiltered edge list dataframe output from `GRNBOOST2` that we used to construct the DEG modules (`MBv_GRNBOOST2_F-test-padj-05_edge-list.csv`) and the TF regulons (`MBv_GRNBOOST2_human-TF_edge-list.csv`).
We also include the results of the `AUCell` estimation of DEG module (`MBv_n119_DEG-modules_AUCell.csv.gz`) or TF regulon (`MBv_n109_TF-regulons_AUCell.csv.gz`) contribution to each spot.
