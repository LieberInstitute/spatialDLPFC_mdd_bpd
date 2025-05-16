setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
	library(scuttle)
	library(ggspavis)
	library(dplyr)
	library(pheatmap)
})
source("code/05_clustering/PRECAST/03-supp_plot-functions.r")
source("code/05_clustering/PRECAST/PRECAST_colorLists.r")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#estimate spatial domain
layer.markers = read.csv("processed-data/04_feature_selection/EXT_TableS9_sig_genes_FDR5perc_enrichment.csv") %>%
  filter(stat>0, spatial_domain_resolution=="Sp09") %>%
  mutate(domain_simple=factor(test, 
                              levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
                              labels=c("Vasc","L1","L2","L3","L4","L5","L6","WM (1)","WM (2)")))

domains = c("Vasc","L1","L2","L3","L4","L5","L6","WM (1)","WM (2)")
### grep'ed for VIM to label L1 (1) as vasc
names(domains) = domains
top100 = lapply(domains, function(x) filter(layer.markers, domain_simple==x) %>% slice_max(n=100, stat) %>% pull(gene))
#pull out only unique markers
t1 = table(unlist(top100))
top100.unique = names(t1)[t1==1]
#length(top100.unique)
#add in WM genes that were in both WM groups
t2 = table(unlist(top100[c("WM (1)","WM (2)")]))
top100.wm = names(t2)[t2==2]
top100.unique = c(top100.unique, top100.wm)
#add GABA genes
gaba = c("SLC6A1","GAD1","GAD2","SST","CCK","TAC1","PVALB","CNR1","RELN")
top100.unique = union(top100.unique, gaba)
cat("\n\nNumber of (mostly) unique layer markers for initial annotation:\n")
length(top100.unique)
#make sure that all genes are present in my dataset
#length(intersect(top100.unique, rowData(spe)$gene_name))==length(top100.unique)
#make df
top100.unique.df = filter(layer.markers, gene %in% top100.unique) %>% group_by(gene) %>% slice_max(n=1, stat) %>% select(gene, domain_simple, fdr, stat) %>%
  mutate(domain_simple=ifelse(gene %in% gaba, "GABA", as.character(domain_simple)))
#add gene_id
lut = rowData(spe)[rowData(spe)$gene_name %in% top100.unique,]
rownames(lut) = lut$gene_name
top100.unique.df$gene_id = lut[top100.unique.df$gene,"gene_id"]

.gene_set = "H-M-markers" #character specifying gene set (starting with 'n' for most instances except for H-M gene set)
.k_clusters = 9 #numeric

lossPlot(.gene_set, .k_clusters)
#if i want to loop this i have to use a for loop so that coldata keeps getting updated
spe <- updateColData(spe, .gene_set, .k_clusters)
quickResaveHDF5SummarizedExperiment(spe)
cat("\nPRECAST clusters", .gene_set, "genes, k=", .k_clusters, "updated to spe with quickResave\n")

uniquepal = precast.colorList[[paste0(.gene_set,"_k",.k_clusters)]][["colors"]]
names(uniquepal) = precast.colorList[[paste0(.gene_set,"_k",.k_clusters)]][["clusters"]]
plist <- generateSpotPlots(spe, .gene_set, .k_clusters, uniquepal)
pdf(file=paste0("plots/05_clustering/PRECAST/PRECAST_", .gene_set, "-k", .k_clusters,"_spot-plots.pdf"), width=12, height=16)
plist[[1]]
plist[[2]]
plist[[3]]
plist[[4]]
plist[[5]]
dev.off()
cat("\nSaved spot plots to:",paste0("plots/05_clustering/PRECAST/PRECAST_", .gene_set, "-k", .k_clusters,"_spot-plots.pdf"), "\n")

annotationHeatmap(spe, .gene_set, .k_clusters, top100.unique.df)

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
