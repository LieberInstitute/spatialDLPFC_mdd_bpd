setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(ggspavis)
	library(ggplot2)
	library(dplyr)
	library(gridExtra)
})

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
#all edges are appropriately identified except for from the following samples
### V13B23-302_A1 and V13B23-302_C1
cdata$true_edges = ifelse(cdata$slide=="V13B23-302", FALSE, cdata$edge_outlier_genes)

#remove areas for low umi
tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>% 
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>% 
  filter(n_spots>5, !is.na(problem_areas_genes.id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$problem_areas_genes.id)
length(unique(remove.areas))


spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Dim spe:",dim(spe),"\n")

stopifnot(identical(rownames(cdata), rownames(colData(spe))))
stopifnot(identical(cdata$in_tissue, spe$in_tissue))
spe$true_edges = cdata$true_edges
spe$problem_areas_genes.id = cdata$problem_areas_genes.id
spe$problem_areas_genes.size = cdata$problem_areas_genes.size


#QC/problem area coldata
spe <- spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")

spe$problem_areas_grouped = ifelse(spe$problem_areas_genes.size<=5, "small", "flag")
spe$problem_areas_grouped = ifelse(spe$problem_areas_genes.size==0, "none", spe$problem_areas_grouped)
colData(spe)[spe$problem_areas_genes.id %in% remove.areas,"problem_areas_grouped"] = "remove"
colData(spe)[spe$true_edges,"problem_areas_grouped"] = "edge"
spe$problem_areas_grouped = as.factor(spe$problem_areas_grouped)
cat("\nspots per problem area group:")
table(spe$problem_areas_grouped)

spe$lowumi = spe$sum_umi<=100
cat("\nlow umi (sum_umi<=100) spots:")
table(spe$lowumi)

#plotting colData
spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

#plotting coldata
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

cat("\nPlotting # genes detected...",format(Sys.time()),"\n")
geneList = lapply(slideList, function(x) {
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="sum_gene", point_size=0.3)+
		scale_color_gradient(low="white", high="navy", labels=function(x) paste0(x/1000,"k"))+
		labs(title=names(l1)[[z]], color="genes")+
		theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
	))
})

cat("\nPlotting chrM ratio...",format(Sys.time()),"\n")
cat("*** Max color limit set to second highest expr_chrM_ratio per sample to help with 100% chrM ratio spots\n")
mitoList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z) {
    max_value = sort(unique(colData(l1[[z]])$expr_chrM_ratio), decreasing=T)[2]
    suppressMessages(plotSpots(l1[[z]], annotate="expr_chrM_ratio", point_size=0.3)+
                       scale_color_gradient(low="white", high="navy", limits=c(0,max_value))+
                       labs(title=names(l1)[[z]], color="chrM")+
                       theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
    )})
})

cat("\nPlotting MBP...",format(Sys.time()),"\n")
mbpList = lapply(slideList, function(x) {
  #cat("Plotting MBP...",format(Sys.time()),"\n")
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="MBP", point_size=0.3, feature_names="gene_name", assay_name="counts")+
                       scale_color_gradient(low="white", high="navy")+
                       labs(title=names(l1)[[z]], color="MBP")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                             panel.background=element_rect(fill="grey30"))
    ))
})

cat("\nPlotting low umi (sum_umi<=100) spots...",format(Sys.time()),"\n")
umiList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="lowumi", point_size=0.3)+
                       scale_color_manual(values=c("grey80","red"))+
                       labs(title=names(l1)[[z]], color="low UMI")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10))
                             #panel.background=element_rect(fill="grey30"), legend.key = element_rect(fill = "white"))
    ))
})

color.palette = c("powderblue","palegoldenrod","#00a000","black","grey80")
names(color.palette) <- c("edge","remove","flag","small","none")

cat("\nPlotting spots by filter group...",format(Sys.time()),"\n")
filtList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z)
    suppressMessages(plotSpots(l1[[z]], annotate="problem_areas_grouped", point_size=0.3)+
                       scale_color_manual(values=color.palette)+
                       labs(title=names(l1)[[z]], color="")+
                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10))
                             #panel.background=element_rect(fill="grey30"), legend.key = element_rect(fill = "white"))
    ))
})

rearrangePlots <- function(slide_id) {
  p0 <- geneList[[slide_id]]
  p1 <- mitoList[[slide_id]]
  p2 <- mbpList[[slide_id]]
  p3 <- umiList[[slide_id]]
  p4 <- filtList[[slide_id]]
  list(p0[[1]], p1[[1]],p2[[1]],p3[[1]],p4[[1]],
       p0[[2]], p1[[2]],p2[[2]],p3[[2]],p4[[2]],
       p0[[3]], p1[[3]],p2[[3]],p3[[3]],p4[[3]],
       p0[[4]], p1[[4]],p2[[4]],p3[[4]],p4[[4]])
}


for(i in names(slideList)) {
  ggsave(file=paste0("plots/03_QC/slide_problem-area_pngs/",i,".png"), 
         do.call(grid.arrange, c(rearrangePlots(i), ncol=5)), 
         bg="white", unit="in", width=15, height=12) 
}

cat("\nSaved to: plots/03_QC/slide_problem-area_pngs/\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()




