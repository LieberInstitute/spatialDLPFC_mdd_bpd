setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(dplyr)
	library(ggplot2)
	library(ggsankey)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")

# load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

# load spatial clustering annotations
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
spotdata = spotdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))
spe$precast_k9_1663 = factor(spotdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI"),
                             labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI"))
spe$smoothed_k9_1663 = factor(spotdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))

# load cell type-driven annotations
res.pc30 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res.pc30 = res.pc30[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(res.pc30)))
spe$seurat_orig = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"),
                          labels=c("Micro.Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"))

spe$seurat_label = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"),
                         labels=c("Micro.Vasc","Astro","L2.3","L2.3","Inhb","L4","L5","L6","Oligo"))

# format for plotting
plot.df = mutate(as.data.frame(colData(spe)), sm_original=paste(precast_k9_1663, "sm og"), smoothed=paste(smoothed_k9_1663, "sm"), 
                 se_original=paste(seurat_orig, "se og"), seurat=paste(seurat_label, "se")) %>%
  make_long(sm_original, smoothed, seurat, se_original)

# prep color palette
cl0 = c(cpList$smoothed.bright, "Vasc"=cpList$transfer.bright[["Micro.Vasc"]], "GABA"=cpList$transfer.bright[["Inhb"]], "low UMI"="grey")
names(cl0) = paste(names(cl0),"sm og")
cl1 = c(cpList$smoothed.bright, "dropped"="grey")
names(cl1) = paste(names(cl1),"sm")
cl2 = cpList$transfer.bright
names(cl2) = paste(names(cl2),"se")
cl3 = cpList$low.res.bright
names(cl3) = paste(names(cl3), "se og")

# minor plotting formatting
plot.df$pretty_node = factor(plot.df$node, levels=c(names(cl0), names(cl1), names(cl2), names(cl3)),
                             labels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA","low UMI",
                                      "L1","L2","L3.4","L5","L6","WM","dropped",
                                      "M.V","Ast","L2.3","L4","L5","L6","Olg","Inb",
                                      "M.V","Ast","L2","L3","L4","L5","L6","Olg","Inb"))

p <- ggplot(plot.df, aes(x=x, next_x=next_x, node=node, next_node=next_node,
                    fill=node, label=pretty_node))+
  geom_sankey()+scale_fill_manual(values=c(cl0, cl1, cl2, cl3), guide="none")+
  geom_sankey_label()+
  theme_minimal()

p_nolabel <- ggplot(plot.df, aes(x=x, next_x=next_x, node=node, next_node=next_node,
                    fill=node, label=pretty_node))+
  geom_sankey()+scale_fill_manual(values=c(cl0, cl1, cl2, cl3), guide="none")+
  #geom_sankey_label()+
  theme_minimal()

pdf(file="plots/publication/supp_clustering_compare/sankey.pdf")
p
p_nolabel
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
