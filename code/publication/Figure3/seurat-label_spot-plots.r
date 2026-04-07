setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(dplyr)
	library(ggplot2)
	library(escheR)
	library(ggrastr)
	library(gridExtra)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")

#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
spotdata = spotdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))
spe$smoothed_k9_1663 = factor(spotdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))

res.pc30 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res.pc30 = res.pc30[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(res.pc30)))
spe$seurat_pc30 = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"),
                           labels=c("Micro.Vasc","Astro","L2.3","L2.3","Inhb","L4","L5","L6","Oligo"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]


pc30.labels = levels(spe$seurat_pc30)
names(pc30.labels) = pc30.labels
plist <- lapply(pc30.labels, function(x) {
  spe_sub$plot.me = spe_sub$seurat_pc30==x
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", point_size=.5) %>%
    add_fill(var="plot.me", point_size=.5)
  p2 <- p+scale_color_manual("", values=cpList$smoothed.bright)+
    scale_fill_manual("", values=c("white","black"))+
    labs(title=x)+theme(text=element_text(size=10), legend.position="none")
  rasterize(p2, layers="Point", dpi=300)
})

ggsave(file="plots/publication/Figure3/seurat-label_spot-plot.pdf", 
       do.call(grid.arrange, c(plist, ncol=4)),
       width=8, height=4)




spe_sub$Inhb = spe_sub$seurat_pc30=="Inhb"
p = make_escheR(spe_sub) %>%
  add_ground(var="smoothed_k9_1663", point_size=.5, stroke=.3) %>%
  add_fill(var="Inhb", point_size=.6)
p1 <- p+scale_color_manual("", values=cpList$smoothed.light)+
  scale_fill_manual("", values=c("white","black"))+
  labs(title="Inhb")+theme(text=element_text(size=10), legend.position="none")

spot.genes = c("CEBPD","APOLD1","SST")
for(i in spot.genes) {
  spe_sub[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

p = make_escheR(spe_sub) %>%
  add_ground(var="smoothed_k9_1663", point_size=.5, stroke=.3) %>%
  add_fill(var="SST", point_size=.6)
p2 <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
 scale_fill_gradient(low="white",high="black", guide="none")+
 labs(title="SST")+
 theme(text=element_text(size=10), plot.title=element_text(face="italic"))

# legend
p3 <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
 scale_fill_gradient(low="white",high="black")+
 theme(text=element_text(size=10), plot.title=element_blank())

pdf(file="plots/publication/Figure3/supp_Inhb-SST_spot-plot.pdf", width=2, height=4)
grid.arrange(rasterize(p1, layers="Point", dpi=300), 
	rasterize(p2, layers="Point", dpi=300), ncol=1)
rasterize(p3, layers="Point", dpi=300)
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
