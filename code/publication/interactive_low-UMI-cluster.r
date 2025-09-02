library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(HDF5Array)
library(escheR)

cpList <- readRDS("plots/colorPalettes.rds")

cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
head(cdata)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
identical(rownames(cdata), colnames(spe))


colData(spe)$smoothed_k9_1663_f = cdata$smoothed_k9_1663_f
spe2 = spe[,!spe$smoothed_k9_1663_f %in% c("Vasc","GABA")]
colData(spe2)$smoothed = factor(spe2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
                               labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))


#standard high quality example section
plist = list(p1, p2)
examples = c("V13B23-329_A1", "V13B23-352_C1", "V13B23-311_C1")


for(i in examples) {
  spe_sub = spe2[,spe2$sample_id==i]
  colData(spe_sub)$is_lowUMI = spe_sub$smoothed=="low UMI"
  colData(spe_sub)$lg10.umi = log10(spe_sub$sum_umi)
  
  p = make_escheR(spe_sub) %>% add_fill(var="smoothed")
  plist[[length(plist)+1]] <- p+scale_fill_manual("PRECAST\n(smoothed)", 
                            values=c(cpList$smoothed.bright,"low UMI"="grey","dropped"="white"))
  
  p = make_escheR(spe_sub) %>% add_fill(var="lg10.umi", point_size = 1.3) %>% 
    add_ground(var="is_lowUMI", point_size=1.3)
  plist[[length(plist)+1]] <- p+scale_color_manual(values=c("grey50","red3"))
}

ggsave(file="plots/publication/examples_lowUMI-cluster_spot-plot.pdf", 
       gridExtra::marrangeGrob(plist, ncol=2, nrow=1, top=NULL),
       width=10, height=5)


#label transfer results
res4 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
identical(rownames(colData(spe)), rownames(res4))

spe$seurat_label = factor(res4$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                            labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))

plist2 = list()
for(i in examples) {
  spe_sub = spe[,spe$sample_id==i]
  colData(spe_sub)$is_lowUMI = spe_sub$smoothed_k9_1663_f=="low UMI"
  colData(spe_sub)$lg10.umi = log10(spe_sub$sum_umi)
  
  p = make_escheR(spe_sub) %>% add_fill(var="lg10.umi", point_size = 1.3) %>% 
    add_ground(var="is_lowUMI", point_size=1.3)
  plist2[[length(plist2)+1]] <- p+scale_color_manual(values=c("grey50","red3"))
  
  
  p = make_escheR(spe_sub) %>% add_fill(var="seurat_label")
  plist2[[length(plist2)+1]] <- p+scale_fill_manual("Label transfer", 
                                                  values=cpList$transfer.bright)
  
}

ggsave(file="plots/publication/examples_lowUMI-cluster_label-transfer_spot-plot.pdf", 
       gridExtra::marrangeGrob(plist2, ncol=2, nrow=1, top=NULL),
       width=10, height=5)
