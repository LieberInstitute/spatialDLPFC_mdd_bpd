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

source("code/05_clustering/PRECAST/04-supp_cluster-names.r")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

for(i in names(annotations)) {
	new.name = paste0(i,"_f")
        colData(spe)[[new.name]] = factor(colData(spe)[[i]], levels= annotations[[i]], labels= names(annotations[[i]]))

}

write.csv(colData(spe)[,c(1,10,12,13,16,26:36)], "processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names = T)
cat("\ncolData csv with cluster results saved to: processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv\n")

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
  #source_dataframe is colData; the rownames need to be the obs names
  #reference type and compare type are columns in source_dataframe
  x_values = unique(source_dataframe[[reference_type]])
  #for all values of reference cluster compute JC
  output_list <- lapply(x_values, function(X) {
    test_x = source_dataframe[[reference_type]]==X
    names(test_x) = rownames(source_dataframe)
    test_y = source_dataframe[[compare_type]]
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
compare.cluster = setdiff(paste0(names(annotations),"_f"), "precast_k9_1663_f")
names(compare.cluster) = compare.cluster
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

plotList = do.call(c, lapply(compare.cluster, function(x) c(plist[x], tilelist[x])))
grobList = marrangeGrob(plotList, nrow=2, ncol=1)
#laymat = rbind(c(1,2,NA),c(3,4,NA))
#pdf(file="plots/05_clustering/PRECAST/compare-clusters_to-k9-1663.pdf", height=6, width=12)
	#grid.arrange(plist[["precast_k5_1663_f"]], plist[["precast_k7_1663_f"]], plist[["precast_k12_1663_f"]], 
        #        tilelist[["precast_k5_1663_f"]], tilelist[["precast_k7_1663_f"]], tilelist[["precast_k12_1663_f"]],
        #        ncol=3)
	#grid.arrange(plist[["precast_k7_HM_f"]], plist[["precast_k9_HM_f"]],
        #        tilelist[["precast_k7_HM_f"]], tilelist[["precast_k9_HM_f"]],
        #        layout_matrix=laymat)
#dev.off()
ggsave("plots/05_clustering/PRECAST/compare-clusters_overlap-with-k9-1663.pdf", grobList, height=6, width=4)
cat("\nOverlap quant. plots saved to: plots/05_clustering/PRECAST/compare-clusters_overlap-with-k9-1663.pdf\n")


umiList <- lapply(names(annotations), function(x) {
	cdata = as.data.frame(colData(spe))
	cdata$plot.me = cdata[[paste0(x,"_f")]]
	ggplot(cdata, aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=sum_umi, fill=plot.me))+
		geom_boxplot(outliers=F)+
		scale_fill_manual(values=annot_colors[[x]])+
		facet_grid(cols=vars(plot.me))+
		scale_y_log10()+labs(x="condition", y="sum UMI (log10 scale)", title=x)+
		theme_minimal()+theme(legend.position="none", panel.grid.major.x=element_blank(),
			text=element_text(size=10), strip.text=element_text(size=10),
                        axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        plot.margin=unit(c(6,12,6,6),"pt"))
})
ggsave("plots/05_clustering/PRECAST/compare-clusters_library-size.png",
	do.call(grid.arrange, c(umiList, ncol=1)),
	bg="white", width=8, height=14)
cat("\nLibrary size boxplots saved to: plots/05_clustering/PRECAST/compare-clusters_library-size.png\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
