setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(ggplot2)
	library(escheR)
	library(ggrastr)
	library(gridExtra)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
#cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]
cdata2 = cdata[cdata$smoothed_k9_1663_f!="low UMI",]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","drop","drop"))

#subset aucell
aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))

aucell = aucell[rownames(colData(spe)),]
m1 = as.matrix(aucell[,1:(grep("seurat_label", colnames(aucell))-1)])

#scale aucell
## new conditional scale functions
find_3MAD = function(distribution) {
  median(distribution)+(3*mad(distribution))
}

find_q99 = function(distribution) {
  quantile(distribution, probs=.995)[[1]]
}

m1 = apply(m1, MARGIN=2, function(x) {
  thresh1 = find_3MAD(x)
  thresh2 = find_q99(x)
  nmax1 = sum(x>thresh1)
  nmax2 = sum(x>thresh2)
  if(nmax2<nmax1) {
    scale_max = thresh2
  } else {
    scale_max = thresh1
  }
  
  # if threshold>max, set threshold to max
  if(scale_max>max(x)) {
    scale_max = max(x)
  }
  
  # normalize to set max
  norm1 = x/scale_max
  norm1[norm1>1] = 1
  return(norm1)
  
})


#add aucell to spe
spot.genes = colnames(m1)
for(i in spot.genes) {
  stopifnot(identical(colnames(spe), rownames(m1)))
  spe[[i]] = m1[,i]
}

  
spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]

mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X","SNHG14","GLUL",
               "CAMK2N1","GAD1","GRIN1","PRKAR1A","UQCRH","APLP1","FTL","PLP1")

cpList <- readRDS("plots/colorPalettes.rds")

plist1 <- lapply(mod_subset, function(x) {
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var=x, point_size = .5)
  p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
    scale_fill_gradient(limits=c(0,1), low="white",high="black", guide="none")+
    labs(subtitle=x)+
    theme(plot.subtitle=element_text(size=9))
  return(rasterize(p2, dpi=200))
})


ggsave(file="plots/publication/Figure3/modules-DEG-subset-refined_NTC-M_spot-plots.pdf",
       arrangeGrob(grobs=plist1, layout_matrix=rbind(1:5,6:10,11:15), top=NULL),
       height=5.5, width=9)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
