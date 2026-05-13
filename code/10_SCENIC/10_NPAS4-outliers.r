setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
  library(gridExtra)
})
set.seed(123)
setAutoBlockSize(1e9)
cpList = readRDS("plots/colorPalettes.rds")

#load aucell
aucell = read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))

#load spe
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




order1 = as.data.frame(colData(spe)) %>% group_by(condition, sex, sample_id, brnum) %>%
  summarise(FOSB_mean=mean(FOSB)) %>%
  ungroup() %>%
  slice_max(n=24, FOSB_mean) %>% arrange(desc(FOSB_mean))

.nrow=6
.ncol=4
lm=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T)


stopifnot(identical(colnames(spe), rownames(aucell)))

for(y in c("FOSB","JUNB","FOS","JUN")) {
  max.valr = 1
  
  #plot the un-normalized version
  #spe[["tmp"]] = aucell[,y]
  #max.val = max(colData(spe)[["tmp"]])
  #max.valr = round(max.val,1)
  #if(max.valr<max.val) max.valr=max.valr+.1
  
  plist1 <- lapply(order1$sample_id, function(x) {
    spe_sub = spe[,spe$sample_id==x]
    p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                           stroke=.3, point_size = .5) |> 
      #add_fill(var="tmp", point_size = .5)
      add_fill(var=y, point_size = .5)
    p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
      scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
      labs(subtitle=unique(spe_sub$sample_id))+
      theme(plot.subtitle=element_text(size=9))
    return(rasterize(p2, dpi=200))
  })
  
  ggsave(file=paste0("plots/10_SCENIC/NPAS4_outliers/top24-FOSB_n119_", y,"_AUCell-regulon-13162-top20_norm_spot-plot.pdf"),
         marrangeGrob(grobs=plist1, layout_matrix=lm, top=paste(y, "(fill scale fixed, can compare across sections)")),
         height=12, width=9)
}


# dlPFC marker genes
spot.genes = c("NPAS4","FOSB","JUNB","EGR4","INHBA","EGR1","DUSP1","BDNF","NR4A3")
for(i in spot.genes) {
  spe[[i]] = logcounts(spe)[rowData(spe)$gene_name==i,]
}

for(y in spot.genes) {
  max.val = max(colData(spe)[[y]])
  max.valr = round(max.val,1)
  if(max.valr<max.val) max.valr=max.valr+.1
  
  plist1 <- lapply(order1$sample_id, function(x) {
    spe_sub = spe[,spe$sample_id==x]
    p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                           stroke=.3, point_size = .5) |> 
      add_fill(var=y, point_size = .5)
    p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
      scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
      labs(subtitle=unique(spe_sub$sample_id))+
      theme(plot.subtitle=element_text(size=9))
    return(rasterize(p2, dpi=200))
  })
  
  ggsave(file=paste0("plots/10_SCENIC/NPAS4_outliers/top24-FOSB_n119_", y, "-expr_AUCell-regulon-13162-top20_spot-plot.pdf"),
         marrangeGrob(grobs=plist1, layout_matrix=lm, top=paste(y, "(fill scale fixed, can compare across sections)")),
         height=12, width=9)
}

#outliers
order1[1:10,c("brnum","sample_id")]

group_by(order1[1:10,], condition, sex) %>% tally()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

