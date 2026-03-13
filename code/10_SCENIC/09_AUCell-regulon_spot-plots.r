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
aucell = read.csv("processed-data/09_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-weighted-regional_AUCell.csv", row.names=1)
cnames = unlist(strsplit(readLines("processed-data/09_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-weighted-regional_AUCell.csv", 1), ","))[-1]
cnames = gsub("\\(\\-\\)", "_enhancer", gsub("\\(\\+\\)", "_promoter", cnames))
colnames(aucell) = cnames

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

#scale aucell
scale_99 = function(distribution) {
  r1 = round(distribution, 4)
  q99 = quantile(r1, probs=c(.99))
  norm1 = r1/q99
  norm1[norm1>1] = 1
  return(norm1)
}
m1 = as.matrix(aucell[,1:(grep("seurat_label", colnames(aucell))-1)])
m1 = apply(m1, MARGIN=2, scale_99)

#add aucell to spe
spot.genes = colnames(m1)
for(i in spot.genes) {
  stopifnot(identical(colnames(spe), rownames(m1)))
  spe[[i]] = m1[,i]
}

#create spe subsets for plotting
#ideal samples
vistoseg.samples = c("308-D1", "329-A1", "342-B1", "332-B1", #NTC M
                     "342-A1", "332-A1", "327-C1", "329-B1",#NTC F
                     "382-C1", "309-D1","352-A1", "329-C1",  #MDD M
                     "382-D1", "352-B1", "309-C1", "380-B1", #"329-D1", #MDD F
                     "382-A1", "352-C1", "342-D1", "327-B1", #BPD M
                     "382-B1", "308-A1", "309-A1", "327-A1") #BPD F

vistoseg.samples = paste0("V13B23-", gsub("-","_", vistoseg.samples))

lowq.samples = c("V13B23-282","V13B23-301","V13B23-302","V13B23-310","V13B23-311","V13B23-403")
lowq.samples = unique(colData(spe)[spe$slide %in% lowq.samples,"sample_id"])

#format for plot layout
.nrow=6
.ncol=4
lm=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T)

for(y in spot.genes) {
  max.valr=1 #for scaled version max is set to 1

  suppressMessages({
  plist1 <- lapply(c(vistoseg.samples, lowq.samples), function(x) {
    spe_sub = spe[,spe$sample_id==x]
    p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                           stroke=.3, point_size = .5) |> 
      add_fill(var=y, point_size = .5)
    p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
      scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
      labs(subtitle=paste(unique(spe_sub$condition), unique(spe_sub$sex)))+
      theme(plot.subtitle=element_text(size=9))
    return(rasterize(p2, dpi=200))
  })
  })

  ggsave(file=paste0("plots/09_SCENIC/AUCell-regulon-13162/example-samples_AUCell-regulon-13162_norm_",gsub("_","-",y),"_spot-plot.pdf"),
         marrangeGrob(grobs=plist1, layout_matrix=lm, top=paste(y, "(fill scale fixed, can compare across sections)")),
         height=12, width=9)
  cat("\nSaved", paste0(y, "...\n"))
}

## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

