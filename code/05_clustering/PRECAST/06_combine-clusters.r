setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(Seurat)
	library(bluster)
	library(ggspavis)
	library(dplyr)
	library(ggplot2)
})

source("code/05_clustering/PRECAST/PRECAST_colorLists.r")

#load in spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

spe$precast_k9_1663_f=factor(spe$precast_k9_1663, levels=precast.colorList[["n1663_k9"]][["clusters"]], labels=precast.colorList[["n1663_k9"]][["annotation"]])
spe$precast_k9_1079_f=factor(spe$precast_k9_1079, levels=precast.colorList[["n1079_k9"]][["clusters"]], labels=precast.colorList[["n1079_k9"]][["annotation"]])

#load in precast object with reduced dimension embeddings
#PRECAST 1663 reduced dims
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1663", withDimnames=F) = reddim[spe$seurat_key,]

# PRECAST 1079 reduced dims
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1079.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1079", withDimnames=F) = reddim[spe$seurat_key,]

# calculate silhouette widths to compare on same graph
set.seed(123)
n1663.results <- as.data.frame(approxSilhouette(reducedDim(spe, "PRECAST_1663"), clusters=spe$precast_k9_1663_f))
set.seed(123)
n1079.results <- as.data.frame(approxSilhouette(reducedDim(spe, "PRECAST_1079"), clusters=spe$precast_k9_1079_f))

comb.df = bind_rows(mutate(tibble::rownames_to_column(n1663.results[,c("cluster","width")]), gene_set="n=1663", cluster=as.character(cluster)),
                    mutate(tibble::rownames_to_column(n1079.results[,c("cluster","width")]), gene_set="n=1079", cluster=ifelse(cluster %in% c("WM 1","WM 2"),"WM",as.character(cluster))
                    )
) %>%
  mutate(cluster=factor(cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM","low UMI")),
         gene_set=factor(gene_set, levels=c("n=1663","n=1079")))

fill_color = precast.colorList[["n1663_k9"]][["colors"]]
names(fill_color) = precast.colorList[["n1663_k9"]][["annotation"]]

p1 <- ggplot(comb.df, aes(x=cluster, y=width, fill=cluster, lty=gene_set))+
  geom_boxplot(outliers=F, position=position_dodge2(preserve="single"))+
  scale_fill_manual(values=fill_color)+
  labs(y="silhouette width", title="Compare PRECAST results: silhouette width")+
  theme_minimal()+theme(legend.position="none", text=element_text(size=10))


# Vasc, L1, and WM all show improved silhouette width with n=1079 results
# look at how expanding these clusters to include n=1079 spots looks spatially
generateSpotPlots <- function(object, test_col, test_colors) {
  #object = spe object to be split into list
  #test_col = colname of colData of variable to plot, will also be root for plot file name
  fname = gsub("\\.","-",test_col)
  #test_colors = group colors for test plot (should be in alphabetical order of groups or named)
  object$dummy_slide = ifelse(object$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", object$slide)
  object$array2 = ifelse(object$slide=="V13B23-283", "A1", object$array)
  seed = levels(as.factor(object$dummy_slide))
  slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
  slideList = lapply(slideList, function(x) {
    do.call(cbind, lapply(x, function(y) object[,object$dummy_slide==y]))
  })
  
  clusPlot = lapply(slideList, function(x) {
    suppressMessages(
      plotSpots(x, annotate=test_col,
                point_size=.3, sample_id="sample_id")+
        scale_color_manual(values=test_colors)+
        facet_grid(rows=vars(dummy_slide), cols=vars(array2))#+
      #theme(panel.background=element_rect(fill="grey30"))
    )
  })
  
  return(clusPlot)
  #pdf(file=paste0("plots/05_clustering/PRECAST_",fname,"_spot-plots.pdf"), width=12, height=16)
  #clusPlot[[1]]
  #clusPlot[[2]]
  #clusPlot[[3]]
  #clusPlot[[4]]
  #clusPlot[[5]]
  #dev.off()
  #cat("\nReassignment spot plots saved to:",paste0("plots/05_clustering/PRECAST_",fname,"_spot-plots.pdf"),"\n")
}


##### combined cluster
## test to see if worth keeping n=1079 L1 as L1 
##### RESOUNDING NO FROM PLOTS
spe$L1.reassign = ifelse(spe$precast_k9_1663_f=="low UMI","low UMI","other")
spe$L1.reassign = ifelse(spe$precast_k9_1663_f=="L1","L1 (n1663)",spe$L1.reassign)
spe$L1.reassign = ifelse(spe$L1.reassign=="other" & spe$precast_k9_1079_f=="L1", "L1 (n1079)", spe$L1.reassign)
table(spe$L1.reassign)
#low UMI L1 (n1079) L1 (n1663)      other 
#  26249      21199      58626     429174
l1_colors = c("#FDBF6F","black","#1F78B4","grey")
names(l1_colors) = c("low UMI","L1 (n1079)","L1 (n1663)","other")
plist <- generateSpotPlots(spe, "L1.reassign", l1_colors)
pdf(file="plots/05_clustering/PRECAST/PRECAST_L1-reassign_spot-plots.pdf", width=12, height=16)
  plist[[1]]
  plist[[2]]
  plist[[3]]
  plist[[4]]
  plist[[5]]
dev.off()
cat("\nReassignment spot plots saved to: plots/05_clustering/PRECAST/PRECAST_L1-reassign_spot-plots.pdf\n")

## test to see if worth keeping n=1079 WM as WM
##### RESOUNDING YES FROM PLOTS
spe$WM.reassign = ifelse(spe$precast_k9_1663_f=="low UMI","low UMI","other")
spe$WM.reassign = ifelse(spe$precast_k9_1663_f=="WM","WM (n1663)",spe$WM.reassign)
spe$WM.reassign = ifelse(spe$precast_k9_1079_f %in% c("WM 1","WM 2") & spe$WM.reassign!="WM (n1663)","WM (n1079)",spe$WM.reassign)
table(spe$WM.reassign) #would increase WM cluster by 1.5
#low UMI      other WM (n1079) WM (n1663) 
#  18570     457305      19338      40035 
wm_colors = c("#FDBF6F","grey","black","#E31A1C")
names(wm_colors) = c("low UMI","other","WM (n1079)","WM (n1663)")
plist <- generateSpotPlots(spe, "WM.reassign", wm_colors)
pdf(file="plots/05_clustering/PRECAST/PRECAST_WM-reassign_spot-plots.pdf", width=12, height=16)
  plist[[1]]
  plist[[2]]
  plist[[3]]
  plist[[4]]
  plist[[5]]
dev.off()
cat("\nReassignment spot plots saved to: plots/05_clustering/PRECAST/PRECAST_WM-reassign_spot-plots.pdf\n")

## test to see if worth keeping n=1079 Vasc as Vasc
##### I SAY YES FROM PLOTS
spe$Vasc.reassign = ifelse(spe$precast_k9_1663_f=="low UMI","low UMI","other")
spe$Vasc.reassign = ifelse(spe$precast_k9_1663_f=="Vasc","Vasc (n1663)",spe$Vasc.reassign)
spe$Vasc.reassign = ifelse(spe$precast_k9_1079_f=="Vasc" & spe$Vasc.reassign!="Vasc (n1663)","Vasc (n1079)",spe$Vasc.reassign)
table(spe$Vasc.reassign) #would double size of vasc cluster
#low UMI        other Vasc (n1079) Vasc (n1663) 
#  24806       484220        14270        11952 
vasc_colors = c("#FDBF6F","grey","black","#FF7F00")
names(vasc_colors) = c("low UMI", "other", "Vasc (n1079)", "Vasc (n1663)")
plist <- generateSpotPlots(spe, "Vasc.reassign", vasc_colors)
pdf(file="plots/05_clustering/PRECAST/PRECAST_Vasc-reassign_spot-plots.pdf", width=12, height=16)
  plist[[1]]
  plist[[2]]
  plist[[3]]
  plist[[4]]
  plist[[5]]
dev.off()
cat("\nReassignment spot plots saved to: plots/05_clustering/PRECAST/PRECAST_Vasc-reassign_spot-plots.pdf\n")


#based on spatial plots and silhouette width, perform semi-supervised expansion/ reassignment/ merging of certain clusters
spe$combined_cluster = as.character(spe$precast_k9_1663_f)
cat("\nStart with PRECAST n=1663 k=9 clusters\n")
(t1= table(spe$combined_cluster))
round(t1/dim(spe)[2], 4)

spe$combined_cluster = ifelse(spe$precast_k9_1079_f %in% c("WM 1","WM 2"), "WM", spe$combined_cluster)
cat("\n\nExpand WM cluster to include all spots in n=1079 WM clusters\n")
(t1 = table(spe$combined_cluster))
round(t1/dim(spe)[2], 4)
#low UMI   GABA     L1     L2     L3     L5     L6   Vasc     WM 
#  18570   6860  58005  96510 139848  70856  73505  11721  59373
spe$combined_cluster = ifelse(spe$precast_k9_1079_f=="Vasc", "Vasc", spe$combined_cluster)
cat("\n\nExpand Vasc cluster to include all spots n=1079 Vasc cluster\n")
(t1 = table(spe$combined_cluster))
round(t1/dim(spe)[2], 4)
#low UMI   GABA     L1     L2     L3     L5     L6   Vasc     WM 
#  17127   6780  50576  95222 137926  70152  72415  25991  59059 

spe_flag = spe[,spe$combined_cluster=="low UMI"]
cat("\n\nIsolate all spots that were part of n=1663 low UMI cluster that have not yet been reassigned\n")
(t1=table(colData(spe_flag)[,c("combined_cluster","precast_k9_1079_f")]))
round(t1/dim(spe)[2], 4)
#..................precast_k9_1079_f
#combined_cluster  Vasc    L1    L2    L3  GABA    L5    L6  WM 1  WM 2
#.........low UMI     0 16313    39   503     0    31   241     0     0

spe_flag2 = spe_flag[,spe_flag$precast_k9_1079_f!="L1"]
dim(spe_flag2)
colData(spe)[rownames(colData(spe_flag2)), "combined_cluster"] = as.character(spe_flag2$precast_k9_1079_f)
cat("\n\nReassign all spots *not* clustered with n=1079 L1 to their n=1079 cluster\n")
(t1 = table(spe$combined_cluster))
round(t1/dim(spe)[2], 4)
#low UMI   GABA     L1     L2     L3     L5     L6   Vasc     WM 
#  16313   6780  50576  95261 138429  70183  72656  25991  59059 

cat("\n\nFinal spot number and proportion for all clusters\n")
spe2 = spe[,spe$combined_cluster!="low UMI"]
(t1 = table(spe2$combined_cluster))
round(t1/dim(spe2)[2], 4)

#add to colData csv
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$combined_cluster = NA
cdata[rownames(colData(spe2)),"combined_cluster"] = spe2$combined_cluster
write.csv(cdata, "processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=T)
cat("\nUpdated PRECAST cluster saved to: processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv\n")

#spot plots of final semi-supervised cluster
spe2$dummy_slide = ifelse(spe2$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe2$slide)
spe2$array2 = ifelse(spe2$slide=="V13B23-283", "A1", spe2$array)
seed = levels(as.factor(spe2$dummy_slide))
slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])

spe2$combined_cluster_f = factor(spe2$combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
slideList = lapply(slideList, function(x) {
  do.call(cbind, lapply(x, function(y) spe2[,spe2$dummy_slide==y]))
})

fill_color = precast.colorList[["n1663_k9"]][["colors"]][1:8]
names(fill_color) = precast.colorList[["n1663_k9"]][["annotation"]][1:8]

clusPlot = lapply(slideList, function(x) {
  suppressMessages(
    plotSpots(x, annotate="combined_cluster_f",
              point_size=.3, sample_id="sample_id")+
      scale_color_manual(values=fill_color)+
      facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
      theme(panel.background=element_rect(fill="grey30"))
  )
})

pdf(file="plots/05_clustering/PRECAST/PRECAST_combined-cluster_spot-plots.pdf", width=12, height=16)
clusPlot[[1]]
clusPlot[[2]]
clusPlot[[3]]
clusPlot[[4]]
clusPlot[[5]]
dev.off()
cat("\nSemi-supervised cluster spot plots saved to: plots/05_clustering/PRECAST/PRECAST_combined-cluster_spot-plots.pdf\n")

#example spot plots
(n16.samples = paste0("V13B23-30",c(1,2,8,9)))
spe_sub2 = spe2[,spe2$slide %in% n16.samples]

# PRECAST 1663 example spot plots
p2 <- plotSpots(spe_sub2, sample_id="sample_id", annotate="combined_cluster_f", point_size=.1)+
  scale_color_manual("combined\ncluster", values=fill_color)+
  facet_wrap(vars(sample_id), ncol=4)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

# combined cluster proportion of spots per sample (tissue composition)
prop.df = mutate(as.data.frame(colData(spe2)), sample_id=as.factor(sample_id)) %>%
  group_by(sample_id, combined_cluster_f, .drop=F) %>% tally() %>% 
  group_by(sample_id) %>% mutate(total_n=sum(n), prop_n=n/total_n) %>%
  #have to left join condition in later otherwise it will keep 119 values for condition with ".drop=F" and that won't do
  left_join(distinct(as.data.frame(colData(spe2)[,c("sample_id","condition")]))) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")))

p3 <- ggplot(prop.df, aes(x=sample_id, y=prop_n, fill=combined_cluster_f))+
  geom_bar(stat="identity", width=.8)+scale_fill_manual(values=fill_color)+
  facet_grid(combined_cluster_f ~ condition, scales="free")+
  theme_minimal()+theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_text(size=6),
                        strip.text.y.right=element_text(angle=0), text=element_text(size=10),
                        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

##sum umi per cluster -- DIDN'T INCLUDE B/C NEXT STEP IS LOOKING AT ALL QC IN PSEUDOBULK
#p4 <- ggplot(as.data.frame(colData(spe2)), aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=sum_umi, 
#           fill=combined_cluster_f))+
#  geom_boxplot(outliers=F)+scale_fill_manual(values=fill_color)+
#  facet_grid(cols=vars(combined_cluster_f))+
#  scale_y_log10()+labs(x="condition")+
#  theme_minimal()+theme(legend.position="none", panel.grid.major.x=element_blank(),
#                        text=element_text(size=10), strip.text=element_text(size=10),
#                        axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
#                        plot.margin=unit(c(6,12,6,6),"pt"))

#look at code in 04_evaluate-clusters_final.r for formatting of above plots
ggsave("plots/05_clustering/PRECAST/PRECAST_combined-cluster_cluster-example-plots.png",
       gridExtra::grid.arrange(p1, p2, p3, layout_matrix=matrix(c(1,2,2,3,3))), bg="white", width=6, height=12, units="in")
cat("\nSaved example plots to: plots/05_clustering/PRECAST/PRECAST_combined-cluster_cluster-example-plots.png\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
