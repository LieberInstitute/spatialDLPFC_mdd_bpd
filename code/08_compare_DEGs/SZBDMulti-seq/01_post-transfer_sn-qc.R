setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

low.res.pal = c("Astro"="#cfa45c","Micro/Vasc"="#911223",
                "Inhb"="#9377AC",
                "L2"="#5D9940","L3"="#5095CD",
                "L4"="#85A0A0",
                "L5"="#ddc94e","L6"="#E45C5F",
                "Oligo"="#D1C4B0")

# load data
cat("\nStart dim BD:\n")
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT.Rdata")
dim(seu_bd)

# save transfer labels
#res1 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-azimuth.csv", row.names=1)
#stopifnot(identical(rownames(res1), colnames(seu_bd)))
res.bd = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(res.bd), colnames(seu_bd)))

seu_bd$seurat_low.res = factor(res.bd$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
seu_bd$prediction.score = res.bd$prediction.score.max

# find clusters to check for single donor contribution
seu_bd <- FindClusters(seu_bd, resolution=0.3)

table(seu_bd@meta.data[,c("individualID","seurat_clusters")])
# three clusters are almost entirely 1 individual (BD4)
# clusters 15 (astro), 8 (L2.3), 14 (olig)

p0 <- ggplot(as.data.frame(seu_bd@meta.data) %>% group_by(seurat_clusters) %>% add_tally(name="n_clus") %>%
	group_by(seurat_clusters, n_clus, individualID) %>% tally() %>% mutate(prop_clus=n/n_clus),
	aes(x=seurat_clusters, y=prop_clus))+
	ggbeeswarm::geom_beeswarm()+scale_y_continuous("prop. cluster single donor", limits=c(0,1))+
	ggtitle("SZBDMulti-seq: BD only")

p1 <- ggplot(as.data.frame(seu_bd@meta.data), aes(x=seurat_clusters, y=prediction.score))+geom_boxplot(outliers=F)+ylim(0,1)#+ggtitle("SZBDMulti-seq: BD only")

p2 <- ggplot(as.data.frame(seu_bd@meta.data), aes(x=seurat_clusters, fill=seurat_low.res))+
        geom_bar(stat="count", position="stack")+
        scale_fill_manual(values=low.res.pal)+theme(legend.position="bottom", legend.key.size=unit(8,"pt"))

# make sure only donors with >1500 nuclei
keep.donors = table(seu_bd$individualID)
keep.donors = names(keep.donors)[keep.donors>1500]
seu_bd = seu_bd[,seu_bd$individualID %in% keep.donors]
# there are no clusters with low prediction scores
dim(seu_bd)

p3 <- ggplot(as.data.frame(seu_bd@meta.data), aes(x=log10(nCount_RNA), group=individualID))+
        stat_ecdf()+facet_wrap(vars(seurat_low.res))+
        geom_hline(aes(yintercept=.5), color="red")+geom_vline(aes(xintercept=3), color="red")

ggsave(file="plots/08_compare_DEGs/SZBDMultiseq_post-sn-qc_bipolar.png",
        gridExtra::grid.arrange(p0, p1, p2, p3, layout_matrix=matrix(c(1,2,3,4,4), ncol=1)),
        width=8, height=11)

#save(seu_bd, file="processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT_post.Rdata")
#cat("\nSeurat object saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT_post.Rdata")


# equivalent control plots
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")

res.c = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-control_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(colnames(seu_con), rownames(res.c)))

seu_con$transfer_low.res = factor(res.c$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
seu_con$prediction.score = res.c$prediction.score.ma

p0 <- ggplot(as.data.frame(seu_con@meta.data) %>% group_by(seurat_clusters) %>% add_tally(name="n_clus") %>%
        group_by(seurat_clusters, n_clus, individualID) %>% tally() %>% mutate(prop_clus=n/n_clus),
        aes(x=seurat_clusters, y=prop_clus))+
        ggbeeswarm::geom_beeswarm()+scale_y_continuous("prop. cluster single donor", limits=c(0,1))+
        ggtitle("SZBDMulti-seq: control only")
p1 <- ggplot(as.data.frame(seu_con@meta.data), aes(x=seurat_clusters, y=prediction.score))+geom_boxplot(outliers=F)+ylim(0,1)
p2 <- ggplot(as.data.frame(seu_con@meta.data), aes(x=seurat_clusters, fill=transfer_low.res))+
        geom_bar(stat="count", position="stack")+
        scale_fill_manual(values=low.res.pal)+theme(legend.position="bottom", legend.key.size=unit(8,"pt"))
p3 <- ggplot(as.data.frame(seu_con@meta.data), aes(x=log10(nCount_RNA), group=individualID))+
        stat_ecdf()+facet_wrap(vars(transfer_low.res))+
        geom_hline(aes(yintercept=.5), color="red")+geom_vline(aes(xintercept=3), color="red")

ggsave(file="plots/08_compare_DEGs/SZBDMultiseq_post-sn-qc_control.png",
        gridExtra::grid.arrange(p0, p1, p2, p3, layout_matrix=matrix(c(1,2,3,4,4), ncol=1)),
        width=8, height=11)

save(seu_con, file="processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_post.Rdata")
cat("\nSeurat object saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_post.Rdata")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
