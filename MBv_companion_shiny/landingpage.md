### Project Overview

We applied spatial transcriptomics (SRT) in the human dorsolateral prefrontal cortex (dlPFC) of 119 donors that were split into diagnostic groups of neurotypical control (NTC), major depressive disorder (MDD), and bipolar disorder (BD). Our donors were selected to be sex- and age-balanced. SRT was achieved with Visium which provides supra-cellular resolution (55um) sequencing linked with spot-specific barcodes that allow for spatial coordinates to be matched with the RNA expression matrix. 

We used two complementary annotation strategies to extract the biological domains of the dlPFC from spot-level gene expression data. Spatial clustering provided domains that correspond to canonical dlPFC cytoarchitectural layers (domain~SP~). Label projection from snRNAseq cell types provided domains that correspond to major cell type groups (domain~CT~) and allow for characterization of trans-laminar cell types like inhibitory neurons. 

We performed pseudobulking for both annotation strategies and used these pseudobulk samples as input for differential expression (DE) testing. Our implementation of two different models allowed us to test for diagnosis-by-sex group differences across the entire dlPFC tissue section (whole-tissue model) and within specific domains (domain-restricted model). The results of these tests are annotated on the violin plots in the "Pseudobulk Expression Plots" tab and browsable in Supplementary Tables 6 and 7 of our manuscript.

We used several analytical tools to gain biological insight from our DE results, including the creation of co-expression modules from DE genes (DEGs). We performed regulon analysis and eQTL analysis and overlapped these results with our DEGs to draw conclusions about how MDD and BD are impacted at the molecular and cellular levels in the dlPFC.

