library(HDF5Array)
library(SpatialExperiment)
library(DelayedArray)
library(escheR)
library(dplyr)
library(ggplot2)

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")



spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
identical(rownames(colData(spe)), rownames(cdata))
spe$precast_k9_1663 = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","L5","L6","WM","GABA","low UMI"),
                             labels=c("Vasc","L1","L2","L3.4","L5","L6","WM","GABA","low UMI"))

spe$smoothed_k9_1663 = factor(cdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))

res.pc30 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res.pc30 = res.pc30[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(res.pc30)))
spe$seurat_pc30 = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                         labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))


spe_sub = spe[,spe$sample_id=="V13B23-329_A1"]

#SVGs
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
avg.expr <- read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1)

#high expr SVGs: AQP4, NEFM, MOBP
#med expr SVGs: HPCAL1, RORB, PCP4
plot.genes = c("GFAP","STMN1","PLP1","HPCAL1","RORB","PCP4")
for(i in plot.genes) {
  colData(spe_sub)[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}


plist <- lapply(plot.genes, function(x) {
  p = make_escheR(spe_sub) %>%
    add_ground(var="in_tissue", stroke=.5, point_size=.8) %>% 
    add_fill(var=x, point_size = .8)
  p+scale_color_manual("", values=c("black"), guide="none")+
    scale_fill_gradient("log2\nCPM", low="white",high="black")+
    labs(title=x)+theme(plot.title=element_text(size=10, face="italic"), 
                        legend.key.size= unit(10, "pt"), legend.title = element_text(size=8),
                        legend.text = element_text(size=7))
})

ggsave(file="plots/publication/example-Br5572_select-SVGs_spot-plots.png", 
       do.call(gridExtra::grid.arrange, c(plist, ncol=3)),
       bg="white", width=6, height=4)



#spatial domains
p = make_escheR(spe_sub) %>% add_fill(var="precast_k9_1663")
p1 <- p+scale_fill_manual("PRECAST\n(original)   ", values=cpList$earthy.pal2)

ggsave(file="plots/publication/example-Br5572_PRECAST-original_spot-plot.png", 
       p1,
       bg="white", width=5, height=5)

p = make_escheR(spe_sub) %>% add_fill(var="smoothed_k9_1663")
p2 <- p+scale_fill_manual("PRECAST\n(smoothed)", values=cpList$smoothed.bright)

ggsave(file="plots/publication/example-Br5572_PRECAST-smoothed_spot-plot.png", 
       p2,
       bg="white", width=5, height=5)


p = make_escheR(spe_sub) %>% add_fill(var="seurat_pc30")
p3 <- p+scale_fill_manual("Seurat (pc30)\nlabel\ntransfer", values=cpList$transfer.bright)

ggsave(file="plots/publication/example-Br5572_seurat-pc30_spot-plot.png", 
       p3,
       bg="white", width=5, height=5)



