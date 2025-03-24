setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scater)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

generate_boxplots <- function(pc.sets, dframe, precast_colors, ylimits=c(-6,6)) {
  #pc.sets is a list of named pc components to plot
  p3.1 <- ggplot(filter(dframe, component %in% pc.sets), aes(x=spot_loading, color=prc_cluster))+
    stat_ecdf()+scale_color_manual(values=precast_colors)+
    facet_wrap(vars(component), ncol=1)+#scale_y_continuous(limits=ylimits, breaks=c(-6,-4,-2,0,2,4,6))+
    coord_cartesian(xlim=ylimits)+
    labs(title="PRECAST cluster", x="zscore of spot loading")+theme_minimal()+
    theme(legend.position="none", plot.title=element_text(hjust=.5),
          panel.grid.minor=element_blank(), #axis.text.x=element_blank(), 
          axis.title.x=element_blank(), #axis.title.x=element_text(size=9), 
          plot.margin = unit(c(6,6,6,0),"pt"), 
          panel.border = element_rect(fill=NA, color="black"))
  p3 <- ggplot(filter(dframe, component %in% pc.sets), aes(x=prc_cluster, y=spot_loading, fill=prc_cluster))+
    geom_boxplot(outliers=F)+scale_fill_manual(values=precast_colors)+
    facet_wrap(vars(component), ncol=1)+#scale_y_continuous(limits=ylimits, breaks=c(-6,-4,-2,0,2,4,6))+
    coord_cartesian(ylim=ylimits)+
    labs(title="PRECAST cluster", y="zscore of spot loading")+theme_minimal()+
    theme(legend.position="none", axis.title.x=element_blank(), plot.title=element_text(hjust=.5),
          panel.grid.minor=element_blank(),
          plot.margin = unit(c(6,6,6,0),"pt"), panel.border = element_rect(fill=NA, color="black"))
  
  p4 <- ggplot(filter(dframe, component %in% pc.sets), aes(x=sample_id, y=spot_loading, color=round))+
    geom_boxplot(outliers=F)+scale_color_manual(values=c("#1b9e77","black","#7570b3"))+
    facet_wrap(vars(component), ncol=1)+#scale_y_continuous(limits=ylimits, breaks=c(-6,-4,-2,0,2,4,6))+
    coord_cartesian(ylim=ylimits)+
    labs(title="sample_id", y="zscore of spot loading", x="color= seq. round")+theme_minimal()+
    theme(legend.position="none", axis.text.x=element_blank(),  plot.title=element_text(hjust=.5),
          panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          axis.title.x=element_text(size=9), plot.margin = unit(c(6,6,6,0),"pt"),
          panel.border = element_rect(fill=NA, color="black"))
  
  p5 <- ggplot(filter(dframe, component %in% pc.sets), aes(x=slide, y=spot_loading, color=round))+
    geom_boxplot(outliers=F)+scale_color_manual(values=c("#1b9e77","black","#7570b3"))+
    facet_wrap(vars(component), ncol=1)+#scale_y_continuous(limits=ylimits, breaks=c(-6,-4,-2,0,2,4,6))+
    coord_cartesian(ylim=ylimits)+
    labs(title="slide", y="zscore of spot loading", x="color= seq. round")+theme_minimal()+
    theme(legend.position="none", axis.text.x=element_blank(), plot.title=element_text(hjust=.5),
          panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          axis.title.x=element_text(size=9), plot.margin = unit(c(6,6,6,0),"pt"),
          panel.border = element_rect(fill=NA, color="black"))
  return(list(prc_ecdf= p3.1, prc_cluster = p3, sample_id = p4, slide = p5))
}

plotPCAs <- function(obj, dimredName, precast_color_list, ylimits=c(-6,6)) {
  #obj = spe object
  #dimred = dimredname
  #precast_color_list is list of 3 vectors: colors, cluster number, and annotation name, with names "colors","clusters","annotation". all 3 need to be in corresponding order and in the order you want plotted on the x axis
  precast_colors = precast_color_list[["colors"]]
  names(precast_colors) = precast_color_list[["annotation"]]
  #set plot titles
  ptitle = unlist(strsplit(dimredName, "_"))
  #get precast var
  prc.var = grep("precast", colnames(colData(obj)), value=T)
  prc.var = grep(ptitle[[length(ptitle)]], prc.var, value=T)
  if(is.null(prc.var)|length(prc.var)!=1) {stop("Trouble finding precast variable in colData")}
  colData(obj)[[prc.var]] = as.factor(colData(obj)[[prc.var]])
  #update ptitle and MNN contingency plan
  if(length(ptitle)==2) {
    if(ptitle[[1]]=="MNN") contingencyPlan=ptitle[[2]] else contingencyPlan=NA
    ptitle = paste0("n=",ptitle[[2]]," ",ptitle[[1]])
  }
  if(length(ptitle)>2) {
    if(ptitle[[1]]=="MNN") contingencyPlan=ptitle[[length(ptitle)]] else contingencyPlan=NA
    ptitle = paste0("n=",ptitle[[length(ptitle)]]," ",ptitle[[1]]," (",paste(ptitle[2:(length(ptitle)-1)], collapse=" "),")")
  }
  #plot percentVar per PC
  if(!is.na(contingencyPlan)) {
    p1 <- ggplot(data.frame("PC"=1:50, "percentVar"=attr(reducedDim(obj, paste0("PCA_",contingencyPlan)),"percentVar")),
                 aes(x=PC, y=percentVar))+
      geom_point(shape=21, size=3, fill="grey")+labs(title=paste0("PCA_",contingencyPlan,": percent Var per ORIGINAL PC"))+theme_bw()+
      theme(panel.grid.minor=element_blank(), text=element_text(size=12))
  } else {
    p1 <- ggplot(data.frame("PC"=1:50, "percentVar"=attr(reducedDim(obj, dimredName),"percentVar")),
                 aes(x=PC, y=percentVar))+
      geom_point(shape=21, size=3, fill="grey")+labs(title=paste0(ptitle,": percent Var per PC"))+theme_bw()+
      theme(panel.grid.minor=element_blank(), text=element_text(size=10))
  }
  #plot percent of variance of PC explained by metadata
  expl.pca = as.data.frame(getExplanatoryPCs(obj, dimredName, n_dimred=50, variables=c("brnum","slide","sum_umi",prc.var))) 
  tmp1 = tidyr::pivot_longer(tibble::rownames_to_column(expl.pca), c("brnum","slide","sum_umi",prc.var), names_to="variable", values_to="var_explained") %>%
    mutate(variable=factor(variable, levels=c("brnum","slide","sum_umi",prc.var)),
           rowname=as.numeric(gsub("PC", "", rowname)))
  p2 <- ggplot(tmp1, aes(x=rowname, y=var_explained, color=variable))+
    geom_point()+geom_line(aes(group=variable))+
    scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","black"))+
    labs(title=paste0(ptitle,": top variables explaining PCs"), x="component", y="% of PC var explained", color="")+
    theme_bw()+
    theme(text=element_text(size=10), legend.position="inside",
          legend.position.inside = c(.8,.5))
  outList = list()
  outList$page1top = list(p1, p2)
  #boxplot of scaled spot loading for 3 major vars
  tmp2 = tidyr::pivot_longer(cbind.data.frame(as.data.frame(colData(obj)[,c("sample_id","slide","round",prc.var)]), scale(reducedDim(obj, dimredName)[,1:25])), 
                                  cols=paste0("PC",1:25), names_to="component", values_to="spot_loading") %>%
    mutate(component=factor(component, levels=paste0("PC",1:25)),
           round=factor(round, levels=c("r1","r2","r2_1","r3"), labels=c("r1","r2","r2","r3")))
  #### to do: make function parameter that is named vector for levels and names so this can be modifiable
  tmp2$prc_cluster = factor(tmp2[[prc.var]], levels=precast_color_list[["clusters"]],
                                 labels=precast_color_list[["annotation"]])
  #generate plots for different sets of PCs
  outList$page1bot = generate_boxplots(paste0("PC",1:5), tmp2, precast_colors, ylimits)
  #outList$page2 = generate_boxplots(paste0("PC",6:15), tmp2, precast_colors, ylimits)
  #outList$page3 = generate_boxplots(paste0("PC",16:25), tmp2, precast_colors, ylimits)
  return(outList)
}


spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
reducedDimNames(spe)

#1663 precast color list
n1663.colorList = list("colors"=c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C","#FDBF6F"),
                       "clusters"=c("5","3","2","7","8","1","6","9","4"),
                       "annotation"=c("Vasc","L1","L2","L3","GABA","L5","L6","WM","flag"))
#1079 precast color list
n1079.colorList = list("colors"=c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C","#8B0000"),
                       "clusters"=c("3","7","5","2","8","9","1","4","6"),
                       "annotation"=c("Vasc","L1","L2","L3","GABA","L5","L6","WM 1","WM 2"))

n1663.plot = plotPCAs(spe, "PCA_1663", n1663.colorList, ylimits=c(-4,4))
ggsave(file="plots/05_clustering/eval-PCA_n1663.png", gridExtra::grid.arrange(n1663.plot[["page1top"]][[1]], n1663.plot[["page1top"]][[2]], 
	n1663.plot[["page1bot"]][[1]], n1663.plot[["page1bot"]][[2]], n1663.plot[["page1bot"]][[3]], n1663.plot[["page1bot"]][[4]],
	layout_matrix=rbind(c(1,1,2,2),c(3,4,5,6),c(3,4,5,6))),
	height=8, width=11, units="in", bg="white")
#do.call(gridExtra::grid.arrange, c(check3[["page2"]], ncol=3))
#do.call(gridExtra::grid.arrange, c(check3[["page3"]], ncol=3))
#dev.off()
cat("\nSaved to: plots/05_clustering/eval-PCA_n1663.png\n")

n1079.plot = plotPCAs(spe, "PCA_1079", n1079.colorList, ylimits=c(-5,5))
ggsave(file="plots/05_clustering/eval-PCA_n1079.png", gridExtra::grid.arrange(n1079.plot[["page1top"]][[1]], n1079.plot[["page1top"]][[2]],
        n1079.plot[["page1bot"]][[1]], n1079.plot[["page1bot"]][[2]], n1079.plot[["page1bot"]][[3]], n1079.plot[["page1bot"]][[4]],
        layout_matrix=rbind(c(1,1,2,2),c(3,4,5,6),c(3,4,5,6))),
        height=8, width=11, units="in", bg="white")
cat("\nSaved to: plots/05_clustering/eval-PCA_n1079.png\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
