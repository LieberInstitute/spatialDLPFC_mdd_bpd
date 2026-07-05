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

# regulons with >10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons1 = filter(regulons, set_size>=10)$TF

#load aucell
aucell= read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)  

# plot avg aucell for each regulon for each sample
df1 = group_by(aucell, condition, sex, sample_id) %>%
  summarise_at(regulons1, mean) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")),
         sex=factor(sex, levels=c("F","M")))
  
df1 = tidyr::pivot_longer(df1, all_of(regulons1), names_to="regulon", values_to="avg_aucell")

ap1.order = filter(ungroup(df1), regulon %in% c("FOS","JUNB","FOSB")) %>% 
  group_by(sample_id) %>% summarise(avg_aucell=mean(avg_aucell)) %>%
  ungroup() %>% mutate(rank_ap1= rank(-avg_aucell)) %>%
  arrange(rank_ap1) %>% pull(sample_id)

top12 = ap1.order[1:12]

#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots
cdata2 = cdata[cdata$smoothed_k9_1663_f!="low UMI",]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","drop","drop"))

spe[["NPAS4"]] = logcounts(spe)[rowData(spe)$gene_name=="NPAS4",]

y="NPAS4"
max.val = max(colData(spe)[[y]])
max.valr = round(max.val,1)
if(max.valr<max.val) max.valr=max.valr+.1


plist1 <- lapply(top12, function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var=y, point_size = .5)
  p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
    scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
    labs(subtitle=unique(spe_sub$brnum))+
    theme(plot.subtitle=element_text(size=9))
  return(rasterize(p2, dpi=200))
})


plist2 <- lapply(top12[1:6], function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var=y, point_size = .5)
  p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
    scale_fill_gradient(limits=c(0,max.valr), low="white",high="black")+
    labs(subtitle=unique(spe_sub$brnum))+
    theme(plot.subtitle=element_text(size=9))
  return(rasterize(p2, dpi=200))
})

pdf(file="plots/publication/supp_regulons-ap1/top-12_ap1-regulon_NPAS4-expr_spot-plot.pdf", height=11, width=4.5)
plot(arrangeGrob(grobs=plist1, ncol=2, top=NULL))
plot(arrangeGrob(grobs=plist2, ncol=1, top="Use for legends"))
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
