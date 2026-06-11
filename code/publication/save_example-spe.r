setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
})

set.seed(123)
setAutoBlockSize(1e9)

cpList = readRDS("plots/colorPalettes.rds")

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#clean up colData (all other demo data available elsewhere)
colData(spe) <- colData(spe)[,c("sample_id","brnum","MBv_sample", #ID keys
                                "condition","sex", #experimental group
                                "array_row","array_col","key","in_tissue", #spatial metrics
                                "sum_umi","sum_gene","expr_chrM_ratio")] #QC metrics

#load domains
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
spe$domain_sp = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","low.UMI","dropped","dropped"))

#add seurat labels
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res = res[rownames(cdata),]
stopifnot(identical(rownames(res), rownames(cdata)))
spe$domain_ct = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                            labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))


#remove low UMI cluster
dim(spe) #28965 535248
spe = spe[,spe$domain_sp!="low.UMI"]
dim(spe) #28965 513200

#extract coldata for GEO
metadata = distinct(as.data.frame(colData(spe))[,c("sample_id","MBv_sample","brnum","condition","sex")])
rownames(metadata) <- NULL
metadata$title = paste(metadata$MBv_sample, gsub("BPD","BD", metadata$condition), metadata$sex, metadata$brnum, sep="_")
metadata$condition <- NULL
metadata$sex <- NULL

write.csv(metadata, "processed-data/publication/GEO_samples_title.csv", row.names=F)

#subset to example samples
example.samples = c("342-A1", "332-A1", "279-A1", "327-C1", #NTC F
                    "382-D1", "023-D1", "334-A1", "309-C1", #MDD F
                    "382-B1", "308-A1", "309-A1", "332-C1", #BPD F
                    "308-D1", "329-A1", "342-B1", "332-B1", #NTC M
                    "382-C1", "309-D1", "340-D1", "328-D1", #MDD M
                    "382-A1", "352-C1", "381-A1", "327-B1") #BPD M

example.samples = paste0("V13B23-", gsub("-","_", example.samples))
example.samples[[6]] = "V13Y10-023_D1"

tmp = spe[,spe$sample_id %in% example.samples]
dim(tmp) #28965 111258
tmp$sample_id = factor(tmp$sample_id, levels=example.samples)
tmp$condition = factor(tmp$condition, levels=c("NTC","MDD","BPD"), labels=c("NTC","MDD","BD"))
tmp$sex = factor(tmp$sex, levels=c("F","M"))

tmp$title = paste(tmp$MBv_sample, tmp$condition, tmp$sex, tmp$brnum, sep="_")

table(distinct(as.data.frame(colData(tmp)[,c("condition","sex","sample_id")]))[,c("condition","sex")])

colData(tmp) <-	colData(tmp)[,c(1:5,15,6:14)]
print(head(colData(tmp)))

#remove counts assay to save space
assays(tmp)$counts <- NULL

regular_matrix_logcounts <- as.matrix(assays(tmp)[["logcounts"]])
sparse_matrix_logcounts <- as(regular_matrix_logcounts, "dgCMatrix")
assays(tmp)$logcounts <- sparse_matrix_logcounts

saveRDS(tmp, "processed-data/publication/spe_n24_example-samples.rds")


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
