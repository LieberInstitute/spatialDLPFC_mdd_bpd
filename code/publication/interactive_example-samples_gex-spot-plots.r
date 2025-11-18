setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
})
set.seed(123)
setAutoBlockSize(1e9)

cpList = readRDS("plots/colorPalettes.rds")

#save name 
save.name = "k9-1663-smooth-final"
fill.palette = cpList$smoothed.bright

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
cat("\nHow many spots are dropped after smoothing clusters:\n")
table(cdata[cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),"smoothed_k9_1663_f"])

cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM"),
	labels=c("L1","L2","L3.4","L5","L6","WM"))


# final spot plots
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

# dlPFC marker genes
spot.genes = c("AQP4","HPCAL1","RORB","PCP4")
## PC3 batch covariate genes
#spot.genes = c("MTRNR2L8","MALAT1","HBA1","MT1X")
for(i in spot.genes) {
  spe[[i]] = logcounts(spe)[rowData(spe)$gene_name==i,]
}

.nrow=6
.ncol=4
lm=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T)

for(y in spot.genes) {
  max.val = max(colData(spe)[[y]])
  max.valr = round(max.val,1)
  if(max.valr<max.val) max.valr=max.valr+.1
  
  plist1 <- lapply(c(vistoseg.samples, lowq.samples), function(x) {
    spe_sub = spe[,spe$sample_id==x]
    p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                           stroke=.3, point_size = .5) |> 
      add_fill(var=y, point_size = .5)
    p2 <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
      scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
      labs(subtitle=paste(unique(spe_sub$condition), unique(spe_sub$sex)))+
      theme(plot.subtitle=element_text(size=9))
    return(rasterize(p2, dpi=200))
  })
  
  ggsave(file=paste0("plots/publication/example-samples_dlPFC-markers_",y,"_spot-plot.pdf"),
	#file=paste0("plots/publication/example-samples_pc3-genes_",y,"_spot-plot.pdf"), 
         marrangeGrob(grobs=plist1, layout_matrix=lm, top=paste(y, "(fill scale fixed, can compare across sections)")),
         height=12, width=9)
}
