# Brief overview of spatial clustering approach

Started with PRECAST due to speed and intrinsic batch correction. Found k=9 worked well (compared to k=7 with several unfinalized feature lists, did not actually compare to k=7 with finalized feature list). Utilizing many feature lists (see `code/04_feature_selection/alt_batch-effect_pipeline`), we consistently found that there is a **low UMI cluster** and we could never identify a L4 cluster. Tried different feature lists and tried k=12. Tried BayesSpace but never worked close to as well as PRECAST. 

### Low UMI cluster
In addition to trying different feature lists, k values, and clustering methods, we also tried two additional methods to address the low UMI cluster (processed in `code/05_clustering/PRECAST/tmp_re-run_PRECAST-n119` scripts; results explored in interactive session).

First, we removed all flagged problem areas and performed PRECAST clustering with n=1663 feature list. This did not improve anything.

Second, we removed all flagged problem areas and also all 3MAD outlier spots that were classified to the low UMI cluster with PRECAST k=9 n1663. We then re-clustered with k=9 n=1663. This did not improve anything.

Inititally we decided to proceed with semi-supervised approach to bolster the Vasc and WM spots that were most impacted by inclusion in the low UMI cluster. 

Upon further reflection (in May 2025), we felt this over-complicated the methods with little gain (e.g. required explanation of BatchSVG, additional filtering of highly expressed SVGs, plus the actualy explanation of the steps to combine the clusters). Since we will perform label transfer anyway, we decided to simplify the methods and just focus on the n1663 feature list results (and also the HM feature list results if that is desired). Similar to the simplification of `code/04_feature_selection`, the plots related to this extended pipeline have been moved to `plots/05_clustering/PRECAST/alt_batch-effect_pipeline`. However, since the evalutation of the different PRECAST results are integrated in the same scripts, we opted not to move these scripts.

## Description of features lists
- *n1663: largest gene set, nnSVG filters of 1) sig in >10 models, 2) sig in >33% of modules run, 3) best rank in any sig model at least 500
- n1629: remove only 34 batch effect genes from n1663
- n1104: remove only top expressing genes with low spcov from n1663
- *n1079: remove both top expressing genes with low spcov AND 34 batch effect genes from n1663
- H-M-markers: Based on DLPFC layer enrichment results in Table S9 from [Huuki-Myers 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11398705/). Top 500 markers for SpD09 results (n=2939 genes).

## PRECAST parameters
- Always 20 iterations (loglik loss plots look good)
- k=9 for all 4 feature lists
- k=12 for main (*) feature lists

## BayesSpace overview and parameters
Many BayesSpace models were tried with different parameters. While there might be able to be more fine-tuning of batch correction, no matter what we were still getting spottier/ less spatially coherent results than with PRECAST. These effects were visible enough that we didn't quantify spatial coherence or contiguity.

I did try using the default BayesSpace parameters (HVGs, BayesSpace-calc PCA, etc.) and that didn't work any better.

### BayesSpace parameters
While not every possible permutation of these parameter combinations were tried, enough models were ran that we got a sense of what was helping (or not). File naming conventions use: 

```
BayesSpace_[reduced dim input]-[feature list]-d[# reduced dims]_q[# of clusters]_iter-[# iterations]_init-[initial clusters]
```

- Reduced dim input: PCA, MNN, PRECAST embeddings tried. Most models used n1079 PCs. See "Batch correction attempts" for more info.
- PCA dims: d=20 (PC1-20), d=13 (PC1-13), and d=12 (PC1-3, PC5-13)
- Iterations: Started with 10K iterations so if `iter-` not specified then 10K was used. For PRECAST we got similar results at 1K, 5K, 10K.
- Initial cluster: Started with PRECAST clusters or kmeans. Either way (and with quite few iterations), there was almost no difference in results despite the initial clusters being wildly different.

### Batch correction attempts

- Evaluation of PCA suggested that there wasn't a strong batch effect so most models focused on PCA to reduce dimensions.
- Fewer n1079 PCs (than n1663 PCs) showed high % of variance explained by brnum/sample_id, slide, and low UMI so these were the default set.
- We did get better performance by optimizing for "high information density" PCs and we looked at whether dropping the single n1079 PC that was highly explained by brnum/sample_id made a difference
- MNN attempted due to the appearance of some sample-specific clusters under some parameters
- Since PRECAST does intrinsic batch correction we tried using the PRECAST embeddings (d=15 used in this case)

