setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(ggspavis)
  library(dplyr)
  library(gridExtra)
  library(scater)
})

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
cat("\nDim spe:",dim(spe),"\n")

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
#all edges are appropriately identified except for from the following samples
### V13B23-302_A1 and V13B23-302_C1
cdata$true_edges = ifelse(cdata$slide=="V13B23-302", FALSE, cdata$edge_outlier_genes)

tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>5, !is.na(problem_areas_genes.id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$problem_areas_genes.id)
length(unique(remove.areas))

stopifnot(identical(rownames(cdata),rownames(colData(spe))))

spe$edge_outlier_genes = cdata$edge_outlier_genes
spe$true_edges = cdata$true_edges

spe$problem_areas_genes.id = cdata$problem_areas_genes.id
spe$problem_areas_genes.size = cdata$problem_areas_genes.size
spe$problem_areas_binary = spe$problem_areas_genes.id %in% remove.areas

spe$lowumi = spe$sum_umi<=100

spe = spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")
cat("Criteria for spot removal\n")
table(colData(spe)[,c("problem_areas_binary","lowumi","true_edges")])

spe$remove_spots = spe$problem_areas_binary | spe$true_edges | spe$lowumi
cat("Total spots for removal:")
table(spe$remove_spots)

spe = spe[,spe$remove_spots==FALSE]
cat("\nDim spe remaining:",dim(spe),"\n")

spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

##QC metrics summary
colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", spe$facet_violin)


colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

seed = levels(as.factor(spe$facet_violin))
slideList = c(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
names(slideList) = slideList
slideList = lapply(slideList, function(x) {
  unlist(lapply(x, function(y)
    sort(unique(colData(spe)[spe$facet_violin==y,"facet_spots"]))
  ))
})

cat("\nPlotting library size...",format(Sys.time()),"\n")
libList = lapply(slideList, function(x) {
  #cat("Plotting library size...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="sum_umi", point_size=0.2)+
                       scale_color_gradient(low="white", high="navy", labels=function(x) paste0(x/1000,"k"))+
                       labs(title=names(l1)[[z]], color="UMI")+
                       theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
    ))
})

cat("\nPlotting # genes detected...",format(Sys.time()),"\n")
geneList = lapply(slideList, function(x) {
  #cat("Plotting genes...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="sum_gene", point_size=0.2)+
                       scale_color_gradient(low="white", high="navy", labels=function(x) paste0(x/1000,"k"))+
                       labs(title=names(l1)[[z]], color="genes")+
                       theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
    ))
})

cat("\nPlotting chrM ratio...",format(Sys.time()),"\n")
cat("*** Max color limit set to second highest expr_chrM_ratio per sample to help with 100% chrM ratio spots\n")
mitoList = lapply(slideList, function(x) {
  #cat("Plotting chrM ratio",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z) {
    max_value = sort(unique(colData(l1[[z]])$expr_chrM_ratio), decreasing=T)[2]
    suppressMessages(plotSpots(l1[[z]], annotate="expr_chrM_ratio", point_size=0.2)+
                       scale_color_gradient(low="white", high="navy", limits=c(0,max_value))+
                       labs(title=names(l1)[[z]], color="chrM")+
                       theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
    )})
})

cat("\nPlotting MBP...",format(Sys.time()),"\n")
markerList1 = lapply(slideList, function(x) {
  #cat("Plotting MBP...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="MBP", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                       scale_color_gradient(low="white", high="navy")+
                       labs(title=names(l1)[[z]], color="MBP")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                             panel.background=element_rect(fill="grey30"))
    ))
})

cat("\nPlotting GAPDH...",format(Sys.time()),"\n")
markerList2 = lapply(slideList, function(x) {
  #cat("Plotting GAPDH...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="GAPDH", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                       scale_color_gradient(low="white", high="navy")+
                       labs(title=names(l1)[[z]], color="GAPDH")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                             panel.background=element_rect(fill="grey30"))
    ))
})

cat("\nPlotting SYT1...",format(Sys.time()),"\n")
markerList3 = lapply(slideList, function(x) {
  #cat("Plotting SYT1...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="SYT1", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                       scale_color_gradient(low="white", high="navy")+
                       labs(title=names(l1)[[z]], color="SYT1")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                             panel.background=element_rect(fill="grey30"))
    ))
})

rearrangePlots <- function(slide_id) {
  p1 <- libList[[slide_id]]
  p2 <- geneList[[slide_id]]
  p3 <- mitoList[[slide_id]]
  p4 <- markerList1[[slide_id]]
  p5 <- markerList2[[slide_id]]
  p6 <- markerList3[[slide_id]]
  list(p1[[1]],p2[[1]],p3[[1]],p4[[1]],p5[[1]],p6[[1]],
       p1[[2]],p2[[2]],p3[[2]],p4[[2]],p5[[2]],p6[[2]],
       p1[[3]],p2[[3]],p3[[3]],p4[[3]],p5[[3]],p6[[3]],
       p1[[4]],p2[[4]],p3[[4]],p4[[4]],p5[[4]],p6[[4]])
}


for(i in names(slideList)) {
  ggsave(file=paste0("plots/03_QC/slide_filtered-QC_pngs/",i,".png"), 
         do.call(grid.arrange, c(rearrangePlots(i), ncol=6)), 
         bg="white", unit="in", width=16, height=12) 
}
cat("\nSaved to: plots/03_QC/slide_filtered-QC_pngs/\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

