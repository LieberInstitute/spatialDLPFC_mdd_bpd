setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(ggspavis)
	library(gridExtra)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
dim(spe)
cdata <- read.csv("processed-data/03_QC/edgeDetection_finalized/colData_found-edges.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#edge-related colData
colData(spe)$umi_3MAD.outlier = cdata$umi_3MAD.outlier
colData(spe)$genes_3MAD.outlier = cdata$genes_3MAD.outlier
colData(spe)$edge_outlier_umi = cdata$edge_outlier_umi
colData(spe)$edge_outlier_genes = cdata$edge_outlier_genes

spe <- spe[,spe$in_tissue]
dim(spe)
colData(spe)$umi_3MAD.outlier <- as.factor(spe$umi_3MAD.outlier)
colData(spe)$genes_3MAD.outlier <- as.factor(spe$genes_3MAD.outlier)
colData(spe)$edge_outlier_umi <- as.factor(spe$edge_outlier_umi)
colData(spe)$edge_outlier_genes <- as.factor(spe$edge_outlier_genes)

#plotting-related colData
spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", colData(spe)$facet_violin)
seed = levels(as.factor(colData(spe)$facet_violin))

colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

#slide/sample list
slideList = c(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
names(slideList) = slideList
slideList = lapply(slideList, function(x) {
        unlist(lapply(x, function(y)
        sort(unique(colData(spe)[spe$facet_violin==y,"facet_spots"]))
        ))
})

#color palettes
color.palette1 = c("grey","red","darkgreen","black")
names(color.palette1) = c("none","slide","sample","both")

color.palette2 = c("grey","red")
names(color.palette2) = c("FALSE","TRUE")

#generate plots
cat("\nPlotting 3MAD UMI outliers...",format(Sys.time()),"\n")
umiList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="umi_3MAD.outlier", point_size=0.2,
              pal=color.palette1)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette1, drop=F)+
      labs(title=names(l1)[[z]], color="UMI")+
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)))
  )
})

cat("\nPlotting 3MAD # genes outliers...",format(Sys.time()),"\n")
genesList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="genes_3MAD.outlier", point_size=0.2,
              pal=color.palette1)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette1, drop=F)+
      labs(title=names(l1)[[z]], color="UMI")+
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)))
  )
})

cat("\nPlotting UMI edges detected...",format(Sys.time()),"\n")
edgeList1 = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="edge_outlier_umi", point_size=0.2, pal=color.palette2)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette2, drop=F)+
      labs(title=names(l1)[[z]], color="UMI")+
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)))
  )
})

cat("\nPlotting # genes edges detected...",format(Sys.time()),"\n")
edgeList2 = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="edge_outlier_genes", point_size=0.2, pal=color.palette2)+
      geom_point(show.legend=TRUE, size=.1)+
      scale_color_manual(values=color.palette2, drop=F)+
      labs(title=names(l1)[[z]], color="genes")+
      theme(legend.text=element_text(size=8), legend.title=element_text(size=10)))
  )
})

cat("Compiling # genes-based edge outlier plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/edgeDetection_finalized/edges-detected_3MAD-genes.pdf", width=12, height=16)
do.call(grid.arrange, c(edgeList2[[1]], ncol=4))
do.call(grid.arrange, c(edgeList2[[2]], ncol=4))
do.call(grid.arrange, c(edgeList2[[3]], ncol=4))
do.call(grid.arrange, c(edgeList2[[4]], ncol=4))
do.call(grid.arrange, c(edgeList2[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/edgeDetection_finalized/edges-detected_3MAD-genes.pdf\n")

rearrangePlots <- function(slide_id) {
  p1 <- umiList[[slide_id]]
  p2 <- genesList[[slide_id]]
  p3 <- edgeList1[[slide_id]]
  p4 <- edgeList2[[slide_id]]
  list(p1[[1]],p2[[1]],p3[[1]],p4[[1]],
       p1[[2]],p2[[2]],p3[[2]],p4[[2]],
       p1[[3]],p2[[3]],p3[[3]],p4[[3]],
       p1[[4]],p2[[4]],p3[[4]],p4[[4]])
}

for(i in names(slideList)) {
        ggsave(file=paste0("plots/03_QC/edgeDetection_finalized/slide_edges_pngs/",i,".png"),
                do.call(grid.arrange, c(rearrangePlots(i), ncol=4)),
                bg="white", unit="in", width=12, height=12)
}
cat("\nSaved to: plots/03_QC/edgeDetection_finalized/slide_edges_pngs/\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
