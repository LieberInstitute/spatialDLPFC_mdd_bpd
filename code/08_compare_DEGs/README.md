### The Goal ###
To quantitatively compare previously published MDD and BD dlPFC datasets with our DE results. To accomplish this we wanted to collect bulk RNA-seq studies to compare with the whole-tissue model and snRNA-seq studies to compare with the domain-restricted model.

Qualifying studies needed to have performed sex-stratified analysis and to have made the full DE results available, including non-significant genes, so that we can correlate logFC or t statistics with our study.

### Bulk RNA-seq ###
For bulk RNA-seq, one BD study fit the bill.
- [Xia et al. 2024](https://doi.org/10.1126/scitranslmed.adh9974) which was a re-analysis of [Gandal et al. 2018](https://doi.org/10.1126/science.aat8127). This study had clear reporting of sample size and sex-specific DE results.

For bulk RNA-seq, two studies included MDD cohorts.
- [Girgenti et al. 2021](https://doi.org/10.1038/s41593-020-00748-7) performed pooled analysis and sex-specific analysis (ST9).
- [Jaffe et al. 2022](https://doi.org/10.1176/appi.ajp.21020162) performed only pooled analysis of M and F together (supplemental CSV, not in excel sheets). 

**Where things stand:** Correlations with our layer-adjusted model, filtered to the 807 F test sig. genes, show good consistency between Girgenti 2021 (MDD) and Xia 2024 (BD).

### snRNA-seq ###
For snRNA-seq, no BD study has been published. However, the psychENCODE release includes a high-quality SZBDMulti-seq dataset (NTC and BD donors) that we have been using for cell type label transfer.
- [Ruzicka et al. 2024](https://doi.org/10.1126/science.adg5136) published on the SCZ portion of this dataset, detailing experimental methods.

For snRNA-seq, two groups have published MDD results, with both groups publishing two studies containing part of the cohort.
- Labonte/Mansouri: The DE results made available are only un-adjusted p-value <.05 and are not cell type-specific. Only half of the full cohort (Labonte-associated data) is available through data repositories.
- Nagy/Maitra: Only adj. p-value <.05 DE results made available, but full cohort is accessible through repositories together under the [Maitra et al. 2023](https://doi.org/10.1038/s41467-023-38530-5) publication. 
- [Chatzinakos et al. 2023](https://doi.org/10.1176/appi.ajp.20220478) did not do sex-stratified analysis but did make all of their DE results available and performed cell type-specific analysis. 

**The proposal:** Use identical Seurat label transfer methods to pseudobulk Maitra 2023 (MDD) and SZBD (BD), and implement DE models in our own hands.

**The nuances:** 
- In addition to removing the mitochondrial genome, the SZBDMulti-seq dataset also removed many highly expressed mitochondrially-destined nuclear genes, including many that we found to be DE and are part of the 1663 SVGs. 
- The SZBD ANOVA for age ~ dx-sex is statistically significant, with BD donors older than NTC donors. Overall the age demographic of this cohort is much older than ours.
- The Maitra MDD dataset did not include donor age (or any continuous donor covariates) as part of sample metadata.
- The Nagy/Maitra MDD samples have the hurdle that the male and female samples were collected and processed at different times for different studies. Many of the technical covariates (batch, chemistry, sequencing platform) are almost impossibly confounded with sex. See the log output of covariate selection step for sex chromosome genes having the most variance explained by technical covariates.

**Where things stand:** 
- Maitra pseudobulking and DE implementation went well and the models looked good (p-value histogram) when including a continuous batch covariate derived from PC2 when 1663 SVGs were used for PCA. Despite identifying a strong DE effect in MDD females and showing most DE results in M/V, Ast, and Olg contexts, these results correlate poorly with our own (different genes).
- The SZBD dataset looks good after pseudobulking, but there is a pervasive hidden batch effect indicated in the p-value histogram that especially effects the males. I'm having a hard time pinning it down, and looked at pseudobulking each donor-cluster by channel as well. Pseudobulk/PCA-wise this doesn't reveal anything obvious, although two PCs (one HVG PCA and one SVG PCA) show a bit better gene-level variance correlation with technical covariates (although introduces its own concerns like pseudobulk samples on <10 cells). However, when I try to run a DE model on this pseudobulk dataset, it takes over 12 hours for voom to run (stopped early), so I'm not sure what is going on.
  - For some reason this is especially bad in L6? I tried looking at the ecdf of prediction score for NTC and BD pre-pseudobulk but I can't see thing obvious.
  - I also explored MetaNeighbor, but the overall consistency is so strong that this didn't show anything helpful. Code in `code/08_compare_DEGs/archive` if you're interested.
  - **Actually the key might be that metadata is wrong. Currently looking at sex chromosome genes!**

**Final thoughts:** One of the models for the Maitra study indicated that some of the tissue damage patterns (MTRNR, HSPA1A, MT1X, MT2A, etc) were significant in the opposite direction. We didn't focus on these results in our study because it seemed quite possible that these DE profiles are driven by peri-mortem factors and a few outlier donors. In our study those outliers were in MDD, it seems very likely that in the Maitra study those outliers are in NTC. If this is an avenue you want to explore, I suggest maybe starting with `BatchSVG`.

