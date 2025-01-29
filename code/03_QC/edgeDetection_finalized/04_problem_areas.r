setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(ggspavis)
	library(raster)
	library(gridExtra)
})
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/edgeDetection_finalized/colData_found-edges.csv", row.names=1)
cdata$umi_3MAD.outlier_binary = cdata$umi_3MAD.outlier_slide | cdata$umi_3MAD.outlier_sample
cdata$umi_3MAD.outlier_binary = ifelse(cdata$in_tissue==FALSE, FALSE, cdata$umi_3MAD.outlier_binary)
cdata$genes_3MAD.outlier_binary = cdata$genes_3MAD.outlier_slide | cdata$genes_3MAD.outlier_sample
cdata$genes_3MAD.outlier_binary = ifelse(cdata$in_tissue==FALSE, FALSE, cdata$genes_3MAD.outlier_binary)
table(cdata[,c("umi_3MAD.outlier_binary","genes_3MAD.outlier_binary")])

sampleList = unique(cdata$sample_id)
names(sampleList) <- sampleList

umi_probs = lapply(sampleList, function(x) {
  problemAreas(cdata[cdata$sample_id==x,c("array_row","array_col","umi_3MAD.outlier_binary")], shifted=F)
})
umi_probs = do.call(rbind, umi_probs)
cdata$problem_areas_umi_size = 0
cdata[umi_probs$spotcode,"problem_areas_umi_size"] = umi_probs$size

genes_probs = lapply(sampleList, function(x) {
  problemAreas(cdata[cdata$sample_id==x,c("array_row","array_col","genes_3MAD.outlier_binary")], shifted=F)
})
genes_probs = do.call(rbind, genes_probs)
cdata$problem_areas_genes_size = 0
cdata[genes_probs$spotcode,"problem_areas_genes_size"] = genes_probs$size

write.csv(cdata, "processed-data/03_QC/edgeDetection_finalized/colData_found-edges_problem-areas.csv", row.names=T)
cat("\ncolData saved to: processed-data/03_QC/edgeDetection_finalized/colData_found-edges_problem-areas.csv\n")


cat("\nPlot # gene-based problem areas...\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
dim(spe)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#problem area-related colData
colData(spe)$problem_areas_genes_size = cdata$problem_areas_genes_size

spe <- spe[,spe$in_tissue]
dim(spe)

#plotting-related colData
spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", colData(spe)$facet_violin)
seed = levels(as.factor(colData(spe)$facet_violin))

colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

slideList2 = list(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList2 = lapply(slideList2, function(x) {
  unlist(lapply(x, function(y)
    sort(unique(colData(spe)[spe$facet_violin==y,"facet_spots"]))
  ))
})

plotList = lapply(slideList2, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="problem_areas_genes_size", point_size=0.2)+
      scale_color_viridis_c(option="F", direction=-1)+
      labs(title=names(l1)[[z]], color="# spots")+
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)))
  )
})

pdf(file="plots/03_QC/edgeDetection_finalized/problem-areas_3MAD-genes.pdf", width=12, height=16)
do.call(grid.arrange, c(plotList[[1]], ncol=4))
do.call(grid.arrange, c(plotList[[2]], ncol=4))
do.call(grid.arrange, c(plotList[[3]], ncol=4))
do.call(grid.arrange, c(plotList[[4]], ncol=4))
do.call(grid.arrange, c(plotList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/edgeDetection_finalized/problem-areas_3MAD-genes.pdf\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
