setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(dplyr)
	library(ggplot2)
	library(escheR)
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
spe$seurat_pc30 = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                           labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))


spe_sub = spe[,spe$sample_id=="V13B23-329_A1"]


pc30.labels = levels(spe$seurat_pc30)
names(pc30.labels) = pc30.labels
plist <- lapply(pc30.labels, function(x) {
  spe_sub$plot.me = spe_sub$seurat_pc30==x
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", point_size=.5) %>%
    add_fill(var="plot.me", point_size=.5)
  p+scale_color_manual("", values=cpList$smoothed.bright)+
    scale_fill_manual("", values=c("white","black"))+
    labs(title=x)+theme(text=element_text(size=10), legend.position="none")
})

ggsave(file="plots/06_pseudobulk/spot-annotation_PRECAST-outline_transfer-fill_spot-plots.png", 
       do.call(grid.arrange, c(plist, ncol=4)),
       bg="white", width=8, height=4)
cat("\nSaved spot-level assignment spot plots to: plots/06_pseudobulk/spot-annotation_PRECAST-outline_transfer-fill_spot-plots.png\n")




spot.genes = c("AQP4","HPCAL1","COL5A2","RORB","PCP4","NR4A2","MOBP","GAD1")
for(i in spot.genes) {
  spe_sub[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

plist2 <- lapply(c("AQP4","HPCAL1","RORB","PCP4"), function(x) {
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", stroke=.5, point_size=.8) %>% 
    add_fill(var=x, point_size = .8)
  p+scale_color_manual("", values=cpList$smoothed.light, guide="none")+
    scale_fill_gradient("log2\nCPM", low="white",high="black")+
    labs(title=x)+theme(plot.title=element_text(size=10, face="italic"), 
                        legend.key.size= unit(10, "pt"), legend.title = element_text(size=8),
                        legend.text = element_text(size=7))
})

ggsave(file="plots/06_pseudobulk/marker-genes_PRECAST-smoothed_spot-plots.png", 
       do.call(grid.arrange, c(plist2, ncol=2)),
       bg="white", width=6, height=6)
cat("\nMarker gene spot plots saved to: plots/06_pseudobulk/marker-genes_PRECAST-smoothed_spot-plots.png\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
