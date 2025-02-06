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

#final cdata
cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=1)
cdata$spotsweeper_outlier = cdata$umi_local.outlier | cdata$genes_local.outlier | cdata$chrM.ratio_local.outlier
cdata[is.na(cdata$spotsweeper_outlier),"spotsweeper_outlier"] = FALSE

cdata$remove_spots = cdata$in_tissue==FALSE | cdata$remove_problem.areas | cdata$spotsweeper_outlier
cdata$problem_area_flag = ifelse(cdata$problem_areas_genes.size>5, TRUE, FALSE)

cat("\n\nRemove spots (off tissue):")
table(cdata[,c("in_tissue","remove_spots")])
tmp = cdata[cdata$in_tissue,]
cat("\n\nRemove in_tissue spots (individual problem area criteria):")
table(tmp[,c("true_edges","remove_spots")])
table(tmp[,c("problem_areas_binary","remove_spots")])
table(tmp[,c("lowumi","remove_spots")])
cat("Remove in_tissue spots (summary problem area criteria):")
table(tmp[,c("remove_problem.areas","remove_spots")])
cat("\n\nRemove in_tissue spots (individual SpotSweeper criteria):")
table(tmp[,c("umi_local.outlier","remove_spots")])
table(tmp[,c("genes_local.outlier","remove_spots")])
table(tmp[,c("chrM.ratio_local.outlier","remove_spots")])
cat("Remove in_tissue spots (summary SpotSweeper criteria):")
table(tmp[,c("spotsweeper_outlier","remove_spots")])

cat("\n\nTotal number of spots removed for any reason:")
table(cdata$remove_spots)
cat("Remaining spots flagged by problem areas:")
table(cdata[cdata$remove_spots==FALSE,"problem_area_flag"])

write.csv(cdata, "processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=T)
cat("\nFinal colData saved to: processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv\n")

#load spe
system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
cat("\nDim spe:",dim(spe),"\n")

stopifnot(identical(rownames(cdata),rownames(colData(spe))))
spe$remove_spots = cdata$remove_spots

#plotting vars
spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", spe$facet_violin)

colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

x_ordered = sort(unique(spe$facet_spots))
colData(spe)$x_facets = ""
colData(spe)[spe$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
colData(spe)[spe$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
colData(spe)[spe$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
colData(spe)[spe$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

#boxplots
df = group_by(as.data.frame(colData(spe)), sample_id) %>% mutate(n_spots_removed=sum(remove_spots)) %>%
  ungroup() %>% filter(remove_spots==FALSE)
cat("Plotting QC boxplots...",format(Sys.time()),"\n")
p1 <- ggplot(df, aes(x=facet_spots, y=sum_umi, fill=n_spots_removed))+
  geom_boxplot(color="grey", outlier.size=.5)+geom_hline(aes(yintercept=1000), lty=2, color="red3")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(100,50000),
                     breaks=c(10^(2:5)), labels=c("100","1k","10k","100k"))+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_c(option="F")+
  labs(x="", y="sum_umi (log10 scale)", title="Library size - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p2 <- ggplot(df, aes(x=facet_spots, y=sum_gene, fill=n_spots_removed))+
  geom_boxplot(color="grey", outlier.size=.5)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_c(option="F")+
  labs(x="", y="sum_gene", title="Detected genes - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p3 <- ggplot(df, aes(x=facet_spots, y=expr_chrM_ratio, fill=n_spots_removed))+
  geom_boxplot(color="grey", outlier.size=.5)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_c(option="F")+
  labs(x="", y="expr_chrM_ratio", title="Mitochondrial fraction - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())

pdf(file="plots/03_QC/filtered_qc-metrics_boxplot.pdf", width=9, height=12)
p1
p2
p3
dev.off()
cat("Saved to: plots/03_QC/filtered_qc-metrics_boxplot.pdf\n")

#by slide spot plots
spe = spe[,spe$remove_spots==FALSE]
cat("\nDim spe remaining:",dim(spe),"\n")

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
mitoList = lapply(slideList, function(x) {
  l1 = x; names(l1) = x
  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
  lapply(1:length(l1), function(z) {
    suppressMessages(plotSpots(l1[[z]], annotate="expr_chrM_ratio", point_size=0.2)+
                       scale_color_gradient(low="white", high="navy")+
                       labs(title=names(l1)[[z]], color="chrM")+
                       theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
    )})
})

#cat("\nPlotting MBP...",format(Sys.time()),"\n")
#markerList1 = lapply(slideList, function(x) {
#  l1 = x; names(l1) = x
#  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
#  
#  lapply(1:length(l1), function(z)
#    suppressMessages(plotSpots(l1[[z]], annotate="MBP", point_size=0.2, feature_names="gene_name", assay_name="counts")+
#                       scale_color_gradient(low="white", high="navy")+
#                       labs(title=names(l1)[[z]], color="MBP")+
#                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
#                             panel.background=element_rect(fill="grey30"))
#    ))
#})

#cat("\nPlotting GAPDH...",format(Sys.time()),"\n")
#markerList2 = lapply(slideList, function(x) {
#  l1 = x; names(l1) = x
#  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
#  
#  lapply(1:length(l1), function(z)
#    suppressMessages(plotSpots(l1[[z]], annotate="GAPDH", point_size=0.2, feature_names="gene_name", assay_name="counts")+
#                       scale_color_gradient(low="white", high="navy")+
#                       labs(title=names(l1)[[z]], color="GAPDH")+
#                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
#                             panel.background=element_rect(fill="grey30"))
#    ))
#})

#cat("\nPlotting SYT1...",format(Sys.time()),"\n")
#markerList3 = lapply(slideList, function(x) {
#  l1 = x; names(l1) = x
#  l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
#  lapply(1:length(l1), function(z)
#    suppressMessages(plotSpots(l1[[z]], annotate="SYT1", point_size=0.2, feature_names="gene_name", assay_name="counts")+
#                       scale_color_gradient(low="white", high="navy")+
#                       labs(title=names(l1)[[z]], color="SYT1")+
#                       theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
#                             panel.background=element_rect(fill="grey30"))
#    ))
#})

rearrangePlots <- function(slide_id) {
  p1 <- libList[[slide_id]]
  p2 <- geneList[[slide_id]]
  p3 <- mitoList[[slide_id]]
  #p4 <- markerList1[[slide_id]]
  #p5 <- markerList2[[slide_id]]
  #p6 <- markerList3[[slide_id]]
  list(p1[[1]],p2[[1]],p3[[1]], #p4[[1]],p5[[1]],p6[[1]],
       p1[[2]],p2[[2]],p3[[2]], #p4[[2]],p5[[2]],p6[[2]],
       p1[[3]],p2[[3]],p3[[3]], #p4[[3]],p5[[3]],p6[[3]],
       p1[[4]],p2[[4]],p3[[4]]) #,p4[[4]],p5[[4]],p6[[4]])
}


for(i in names(slideList)) {
  ggsave(file=paste0("plots/03_QC/slide_filtered-QC_pngs/",i,".png"), 
         do.call(grid.arrange, c(rearrangePlots(i), ncol=3)), 
         bg="white", unit="in", width=8, height=12) 
}
cat("\nSaved to: plots/03_QC/slide_filtered-QC_pngs/\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

