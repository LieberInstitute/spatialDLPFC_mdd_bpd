setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(edgeR)
	library(scater)
})

# the normalized pseudobulk version is already filtered to exclude some of the top marker genes so load the pre-normalized spe and normalize
#load("processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-with-lowUMI_sample-n1663-k9_norm.Rdata")
#length(intersect(top.genes, rownames(spe_pseudo))) #728
load("processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-with-lowUMI_sample-n1663-k9.Rdata")
#length(intersect(top.genes, rownames(spe_pseudo))) #752

tmp = calcNormFactors(spe_pseudo)
x = cpm(tmp, log=T, prior.count=2)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x

## reverse the factor order for plotting
#spe_pseudo$precast_k9_1663_rev = factor(spe_pseudo$precast_k9_1663, levels=rev(c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI")),
#                             labels=rev(c("Mng","L1","L2","L3.4","GABA","L5","L6","WM","low UMI")))

# load in spatialDLPFC enrichment results
layer_modeling_results <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_modeling_results")
t1 = layer_modeling_results$enrichment[,c(grep("t_stat",colnames(layer_modeling_results$enrichment), value=T),"ensembl")]
fix.names = c("Sp09D01"="Mng","Sp09D02"="L1","Sp09D03"="L2","Sp09D04"="L5","Sp09D05"="L3","Sp09D06"="WM.1","Sp09D07"="L6","Sp09D08"="L4","Sp09D09"="WM.2")
m1 = as.matrix(t1[,1:9])
colnames(m1) <- fix.names
rownames(m1) = t1$ensembl
#dim(m1)

# collect top 100 from each domain
source("code/06_pseudobulk/custom_functions.r")
top.genes = getTopGenes(m1, top_n=100)
length(top.genes) #752

top.mtx = getTopGenes(m1, top_n=100, return_matrix=T)
top.mtx = top.mtx[top.genes,]

# supplement with GABA
top.mtx = cbind(top.mtx, "GABA"=rep(FALSE, length.out=nrow(top.mtx)))

inhb.genes = c("SLC6A1","GAD1","GAD2","SST","CCK","CNR1","TAC1","PVALB","RELN","SLC32A1")
inhb.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% inhb.genes]
top.mtx = rbind(top.mtx, 
      matrix(c(rep(FALSE, length.out=ncol(top.mtx)-1), TRUE), ncol=ncol(top.mtx), nrow=length(setdiff(inhb.ids, top.genes)),
       byrow=T, dimnames = list(setdiff(inhb.ids, top.genes), colnames(top.mtx)))
)

# formatting for heatmap
row_annot = as.data.frame(apply(top.mtx, 2, as.character))
rownames(row_annot) = rownames(top.mtx)

cpList <- readRDS("plots/colorPalettes.rds")
annot_colors = lapply(colnames(row_annot), function(x) {
  x1 = factor(x, levels=c("Mng","L1","L2","L3","L4","GABA","L5","L6","WM.1","WM.2"),
              labels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo","Oligo"))
  
  c("FALSE"="white","TRUE"=cpList$low.res.bright[[as.character(x1)]])
  })
names(annot_colors) = colnames(row_annot)

# make heatmap
phm = plotGroupedHeatmap(spe_pseudo, features=rownames(row_annot), group="precast_k9_1663", 
                           color=colorRampPalette(c("white","grey90","grey70","black"))(100),
                           scale=T, center=T, show_rownames=F, cluster_col=F, angle_col=90,
                           annotation_row=row_annot[,rev(c("Mng","L1","L2","L3","L4","GABA","L5","L6","WM.1","WM.2"))], 
                           annotation_colors=annot_colors)

pdf(file="plots/publication/supp_clustering_PRECAST/spatialDLPFC-markers_initial-annotation-heatmap.pdf")
plot(phm[[4]])
dev.off()



cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
