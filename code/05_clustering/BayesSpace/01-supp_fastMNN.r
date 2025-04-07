#################################################
####### this was run on my local computer while JHPCE was down so there is no .sh file'
#################################################

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(batchelor)
  library(BiocParallel)
  library(BiocSingular)
  library(dplyr)
})
set.seed(123)

#spe <- loadHDF5SummarizedExperiment(dir="MBv_local/", prefix="spe_n120_postQC_norm_")
load("MBv_local/new_files/spe_PCA_MNN_empty-assays.Rdata")
reducedDimNames(dummy_spe)

#old.geneList <- readRDS("MBv_local/nnSVG-eval_geneList.rds")
#top.1054 = setdiff(old.geneList$qual_genes, old.geneList$top.decile_low.spcov)
#id.1054 = rownames(spe)[rowData(spe)$gene_name %in% top.1054]
#length(id.1054)

#specify merge order
samp.data = distinct(as.data.frame(colData(dummy_spe)[,c("sample_id","brnum","condition","sex","round","age","PMI","RIN")]))
ntcM = filter(samp.data, condition=="NTC", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
ntcF = filter(samp.data, condition=="NTC", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
mddM = filter(samp.data, condition=="MDD", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
mddF = filter(samp.data, condition=="MDD", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r2_1","r1","r3"))) %>% arrange(round) %>% pull(brnum)
bpdM = filter(samp.data, condition=="BPD", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
bpdF = filter(samp.data, condition=="BPD", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
m.order = list(ntcM, ntcF, mddM, mddF, bpdM, bpdF)

cat("\nStart fastMNN...\n")
Sys.time()
#mnn <- fastMNN(spe, batch=spe$brnum, cos.norm = T, get.variance = T, subset.row= id.1054, BSPARAM=RandomParam())
#mnn <- reducedMNN(reducedDim(dummy_spe,'PCA_1079'), batch=dummy_spe$brnum, k=10, merge.order=m.order,
#                  BPPARAM=SerialParam())
mnn <- reducedMNN(reducedDim(dummy_spe,'PCA_1079'), batch=dummy_spe$slide, k=10, #merge.order=m.order,
                  BPPARAM=SerialParam())
#reducedDim(dummy_spe,'MNN_1638') <- mnn$corrected
#reducedDimNames(dummy_spe)
Sys.time()
#save(dummy_spe, file="MBv_local/new_files/spe_PCA_MNN_empty-assays.Rdata")
#cat("\nSaved to:MBv_local/new_files/spe_PCA_MNN_empty-assays.Rdata")
save(mnn, file="MBv_local/new_files/spe_fastMNN_n1079-k10_slide.Rdata")
cat("Saved to: MBv_local/new_files/spe_fastMNN_n1079-k10_slide.Rdata")

