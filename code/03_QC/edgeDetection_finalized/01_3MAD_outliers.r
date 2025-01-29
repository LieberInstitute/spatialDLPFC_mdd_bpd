setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(ggspavis)
	library(scuttle)
})

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
dim(spe)

spe$lg10.umi = log10(spe$sum_umi)
spe$lg10.genes = log10(spe$sum_gene)

cat("\n3MAD outlier to ID poor quality edges:", format(Sys.time()),"\n")
colData(spe)$umi_3MAD.outlier_slide = isOutlier(spe$lg10.umi, subset=spe$in_tissue, batch=spe$slide, type="lower", nmads=3)
colData(spe)$genes_3MAD.outlier_slide = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$slide, type="lower", nmads=3)

colData(spe)$umi_3MAD.outlier_sample = isOutlier(spe$lg10.umi, subset=spe$in_tissue, batch=spe$sample_id, type="lower", nmads=3)
colData(spe)$genes_3MAD.outlier_sample = isOutlier(spe$lg10.genes, subset=spe$in_tissue, batch=spe$sample_id, type="lower", nmads=3)

colData(spe)$umi_3MAD.outlier = ifelse(spe$in_tissue==FALSE, "off tissue", paste(spe$umi_3MAD.outlier_sample, spe$umi_3MAD.outlier_slide))
colData(spe)$umi_3MAD.outlier = factor(spe$umi_3MAD.outlier, 
	levels=c("off tissue","FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("off tissue","none","slide","sample","both"))

colData(spe)$genes_3MAD.outlier = ifelse(spe$in_tissue==FALSE, "off tissue", paste(spe$genes_3MAD.outlier_sample, spe$genes_3MAD.outlier_slide))
colData(spe)$genes_3MAD.outlier = factor(spe$genes_3MAD.outlier, 
        levels=c("off tissue","FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
        labels=c("off tissue","none","slide","sample","both"))

write.csv(colData(spe), "processed-data/03_QC/edgeDetection_finalized/colData_3MAD.csv", row.names=T)

#plot
colData(spe)$facet_spots = paste(spe$round, spe$sample_id)

spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

spe$seed = paste(spe$round, spe$dummy_slide)
spe$seed = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", spe$seed)
seed = levels(as.factor(spe$seed))
slideList = list(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList = lapply(slideList, function(x)
  unlist(lapply(x, function(y)
    sort(unique(colData(spe)[spe$seed==y,"facet_spots"]))
  ))
)

color.palette = c("grey","white","red","darkgreen","black")
names(color.palette) = c("off tissue","none","slide","sample","both")

plotList = lapply(slideList, function(x) {
  cat("Plotting...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="umi_3MAD.outlier", in_tissue=NULL, point_size=0.2,
              pal=color.palette)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette, drop=F)+
      labs(title=names(l1)[[z]], color="UMI outlier"))
  )
})

cat("Compiling UMI outlier plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/edgeDetection_finalized/3MAD-outliers_umi.pdf", width=12, height=16)
        PRECAST::drawFigs(plotList[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
cat("\nSaved here: plots/03_QC/edgeDetection_finalized/3MAD-outliers_umi.pdf\n")

plotList = lapply(slideList, function(x) {
  cat("Plotting...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="genes_3MAD.outlier", in_tissue=NULL, point_size=0.2,
              pal=color.palette)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette, drop=F)+
      labs(title=names(l1)[[z]], color="# genes outlier"))
  )
})

cat("Compiling # genes outlier plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/edgeDetection_finalized/3MAD-outliers_genes.pdf", width=12, height=16)
        PRECAST::drawFigs(plotList[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
        PRECAST::drawFigs(plotList[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
cat("\nSaved here: plots/03_QC/edgeDetection_finalized/3MAD-outliers_genes.pdf\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
