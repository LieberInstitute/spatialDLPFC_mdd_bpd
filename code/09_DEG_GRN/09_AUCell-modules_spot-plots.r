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

REGULON_TYPE = "modules-DEG-subset-refined"
REG_KEY = "DEG-subset-refined"

#load aucell
aucell = read.csv(paste0("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_", REGULON_TYPE,"_AUCell.csv"), row.names=1)
colnames(aucell) = gsub("Regulon\\.for\\.", "", colnames(aucell))

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


plot_expr = TRUE
#plot_expr = FALSE

if(!plot_expr) {
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

} else {
	spot.genes = colnames(m1)
	for(i in spot.genes) {
		spe[[i]] = logcounts(spe)[rowData(spe)$gene_name==i,]
	}

}



#create spe subsets for plotting
#ideal samples
#subset to example samples
example.samples = c("342-A1", "332-A1", "279-A1", "327-C1", #NTC F
                    "382-D1", "023-D1", "334-A1", "309-C1", #MDD F
                    "382-B1", "308-A1", "309-A1", "332-C1", #BPD F
                    "308-D1", "329-A1", "342-B1", "332-B1", #NTC M
                    "382-C1", "309-D1", "340-D1", "328-D1", #MDD M
                    "382-A1", "352-C1", "381-A1", "327-B1") #BPD M

example.samples = paste0("V13B23-", gsub("-","_", example.samples))
example.samples[[6]] = "V13Y10-023_D1"

lowq.samples = c("V13B23-282","V13B23-301","V13B23-302","V13B23-310","V13B23-311","V13B23-403")
lowq.samples = unique(colData(spe)[spe$slide %in% lowq.samples,"sample_id"])

#format for plot layout
.nrow=6
.ncol=4
lm=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T)

for(y in spot.genes) {
  if(!plot_expr) {
	max.valr=1 #for scaled version max is set to 1
  } else {
	max.val = max(colData(spe)[[y]])
	max.valr = round(max.val,1)
	if(max.valr<max.val) max.valr=max.valr+.1
  }
  
  suppressMessages({
  plist1 <- lapply(c(example.samples, lowq.samples), function(x) {
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

  if(!plot_expr) file_name = paste0("plots/09_DEG_GRN/AUCell-modules/example-samples_",gsub("_","-",y),"_AUCell-modules-13162-", REG_KEY, "_spot-plot.pdf")
  if(plot_expr) file_name = paste0("plots/09_DEG_GRN/AUCell-modules/example-samples_",gsub("_","-",y),"-expr_AUCell-modules-13162-", REG_KEY, "_spot-plot.pdf")

  ggsave(file=file_name,
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

