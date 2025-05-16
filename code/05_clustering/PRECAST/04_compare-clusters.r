setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(mclust)
	library(ggplot2)
	library(dplyr)
	library(gridExtra)
})

source("code/05_clustering/PRECAST/PRECAST_colorLists.r")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

spe$precast_k9_1663_f=factor(spe$precast_k9_1663, levels=precast.colorList[["n1663_k9"]][["clusters"]], labels=precast.colorList[["n1663_k9"]][["annotation"]])
spe$precast_k9_1079_f=factor(spe$precast_k9_1079, levels=precast.colorList[["n1079_k9"]][["clusters"]], labels=precast.colorList[["n1079_k9"]][["annotation"]])
spe$precast_k9_1629_f=factor(spe$precast_k9_1629, levels=precast.colorList[["n1629_k9"]][["clusters"]], labels=precast.colorList[["n1629_k9"]][["annotation"]])
spe$precast_k9_1104_f=factor(spe$precast_k9_1104, levels=precast.colorList[["n1104_k9"]][["clusters"]], labels=precast.colorList[["n1104_k9"]][["annotation"]])
spe$precast_k12_1079_f=factor(spe$precast_k12_1079, levels=precast.colorList[["n1079_k12"]][["clusters"]], labels=precast.colorList[["n1079_k12"]][["annotation"]])
spe$precast_k12_1663_f=factor(spe$precast_k12_1663, levels=precast.colorList[["n1663_k12"]][["clusters"]], labels=precast.colorList[["n1663_k12"]][["annotation"]])
spe$precast_k9_HM_f=factor(spe$precast_k9_HM, levels=precast.colorList[["H-M-markers_k9"]][["clusters"]], labels=precast.colorList[["H-M-markers_k9"]][["annotation"]])

#write.csv(colData(spe)[,c("sample_id","slide","array","brnum","age","PMI","RIN",
#	"precast_k9_1663","precast_k9_1663_f","precast_k9_1079","precast_k9_1079_f",
#	"precast_k9_1629","precast_k9_1629_f","precast_k9_1104","precast_k9_1104_f",
#	"precast_k12_1663","precast_k12_1079","precast_k9_HM")], 
#	"processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names = T)
cat("\ncolData csv with cluster results saved to: processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv\n")

#keep in code just for reference
#adjustedRandIndex(spe$precast_k9_1663, spe$precast_k9_1629)
##0.9206876 #remove just batch genes
#adjustedRandIndex(spe$precast_k9_1663, spe$precast_k9_1104)
##0.4861294 #removing just high expr low spcov
#adjustedRandIndex(spe$precast_k9_1663, spe$precast_k9_1079)
##0.6955969 #removing both high expr low spcov AND batch genes

#adjustedRandIndex(spe$precast_k9_1079, spe$precast_k9_1629)
##0.7029754 #adding just high expr low spcov
#adjustedRandIndex(spe$precast_k9_1079, spe$precast_k9_1104)
##0.589903 #adding just batch genes

################# pairwise jaccard
jcoef <- function(x, y) { 
  # x is a named T/F vector for reference cluster ID, with names being spot codes
  # y is a named factor, numeric, or character vector of all comparison cluster IDs, with names being spotcodes
  if(class(y)!="factor") y=as.factor(y)
  x.ids = names(x)[x]
  y.list = levels(y)
  names(y.list) = y.list
  #for all levels of comparison factor, return JC with reference
  sapply(y.list, function(z) {
    z.ids = names(y)[y==z]
    length(intersect(x.ids, z.ids))/length(union(x.ids, z.ids))
  })
}

pairwise_jc <- function(source_dataframe, reference_type, compare_type) {
  #source_dataframe should be pretty self_explanatory, but the rownames need to be the obs names
  #reference type and compare type are columns in match_manual
  x_values = unique(source_dataframe[[reference_type]])
  #for all values of reference cluster compute JC
  output_list <- lapply(x_values, function(X) {
    test_x = source_dataframe[[reference_type]]==X
    #names(test_x) = source_dataframe$spot_id
    names(test_x) = rownames(source_dataframe)
    test_y = source_dataframe[[compare_type]]
    #names(test_y) = source_dataframe$spot_id
    names(test_y) = rownames(source_dataframe)
    #for reference cluster ID X, compute JC for all comparison clusters
    jc_output <- jcoef(test_x, test_y)
    cbind.data.frame(ref_type=rep(reference_type, length(jc_output)),
                     ref_clus=rep(X, length(jc_output)),
                     comp_type=rep(compare_type, length(jc_output)),
                     comp_clus=names(jc_output),
                     j.coef=as.numeric(jc_output))
  })
  do.call(rbind, output_list)
}

#generate plots 
compare.cluster = c("precast_k9_1629_f", "precast_k9_1104_f", "precast_k9_1079_f","precast_k12_1663_f","precast_k12_1079_f","precast_k9_HM_f")
#jaccard plots
plist <- lapply(compare.cluster, function(X) {
  ggplot(pairwise_jc(colData(spe), "precast_k9_1663_f", X), 
         aes(x=ref_clus, y=j.coef, label=as.character(comp_clus)))+
    geom_text(size=3)+ylim(0,1)+scale_x_discrete(expand=expansion(mult=.1))+
    labs(title=paste0("vs. ", X, " (ARI= ", round(adjustedRandIndex(spe$precast_k9_1663, colData(spe)[[X]]), 3),")"), 
         x="PRECAST n=1663 k=9", y="Jaccard coef.")+
    theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=8))
})
#heatmap plots
tilelist <- lapply(compare.cluster, function(X) {
  tmp.df = as.data.frame(table(colData(spe)[,c("precast_k9_1663_f",X)]))
  colnames(tmp.df) = c("reference","query","Freq")
  tmp.df = group_by(tmp.df, reference) %>% mutate(Total=sum(Freq)) %>% ungroup() %>%
    mutate(Prop=Freq/Total, reference=factor(reference, levels=levels(spe$precast_k9_1663_f)),
           query=factor(query, levels=rev(levels(colData(spe)[[X]]))))
  ggplot(tmp.df, aes(y=query, x=reference, fill=Prop))+
    geom_tile(color="grey50", linewidth=.3)+scale_fill_gradient(low="white",high="black", limits=c(0,1))+
    geom_text(data=union(group_by(tmp.df, query) %>% slice_max(n=2, Freq), group_by(tmp.df, reference) %>% slice_max(n=2, Freq)) %>% 
                filter(Freq>1000, Prop>.1) %>%
                mutate(text_value= paste0(round(Freq/1000, 1), "k")), 
              aes(label=text_value), color="red", size=2, fontface="bold")+
    labs(x="PRECAST n=1663 k=9", y=X, fill="Prop. of\nn1663\ncluster")+
    theme_minimal()+theme(text=element_text(size=8), legend.key.size=unit(6,"pt"),
		legend.title=element_text(margin=margin(0,0,4,0,"pt")),
		legend.box.spacing = unit(2,"pt"),
		legend.margin=margin(0,0,0,0,"pt"),
		legend.box.margin = margin(0,2,0,2,"pt"))
})


ggsave("plots/05_clustering/PRECAST/compare-clusters_1663-1079-k9-k12.png",
	grid.arrange(plist[[3]], plist[[4]], plist[[5]], 
		tilelist[[3]], tilelist[[4]], tilelist[[5]],
		ncol=3),
	bg="white", height=6, width=12, units="in"
)
cat("\nSaved to: plots/05_clustering/PRECAST/compare-clusters_1663-1079-k9-k12.png\n")

ggsave("plots/05_clustering/PRECAST/compare-clusters_1629-1104-HM-k9.png",
        grid.arrange(plist[[1]], plist[[2]], plist[[6]],
                tilelist[[1]], tilelist[[2]], tilelist[[6]],
                ncol=3),
        bg="white", height=6, width=12, units="in" 
)
cat("\nSaved to: plots/05_clustering/PRECAST/compare-clusters_1629-1104-HM-k9.png\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
