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

load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT.Rdata")
cat("\nStart dim males:\n")
dim(seu_m)
load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT.Rdata")
cat("\nStart dim females:\n")
dim(seu_f)

res.m = read.csv("processed-data/08_compare_DEGs/Maitra/results_label-transfer_MBv-filtered_ref-SZBD-control_query-male_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
res.f = read.csv("processed-data/08_compare_DEGs/Maitra/results_label-transfer_MBv-filtered_ref-SZBD-control_query-female_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)

stopifnot(identical(colnames(seu_m), rownames(res.m)))
stopifnot(identical(colnames(seu_f), rownames(res.f)))

seu_m$seurat_low.res = factor(res.m$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
seu_m$prediction.score = res.m$prediction.score.max
seu_f$seurat_low.res = factor(res.f$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
seu_f$prediction.score = res.f$prediction.score.max

# individual donors per cluster
p0.m <- ggplot(as.data.frame(seu_m@meta.data) %>% group_by(seurat_clusters) %>% add_tally(name="n_clus") %>%
        group_by(seurat_clusters, n_clus, Sample) %>% tally() %>% mutate(prop_clus=n/n_clus),
        aes(x=seurat_clusters, y=prop_clus))+
        ggbeeswarm::geom_beeswarm()+scale_y_continuous("prop. cluster single donor", limits=c(0,1))+
        ggtitle("Maitra 2023: males")
p0.f <- ggplot(as.data.frame(seu_f@meta.data) %>% group_by(seurat_clusters) %>% add_tally(name="n_clus") %>%
        group_by(seurat_clusters, n_clus, Sample) %>% tally() %>% mutate(prop_clus=n/n_clus),
        aes(x=seurat_clusters, y=prop_clus))+
        ggbeeswarm::geom_beeswarm()+scale_y_continuous("prop. cluster single donor", limits=c(0,1))+
        ggtitle("Maitra 2023: females")

p1.m <- ggplot(as.data.frame(seu_m@meta.data), aes(x=seurat_clusters, y=prediction.score))+geom_boxplot(outliers=F)+ylim(0,1)#+ggtitle("Maitra 2023: males")
p1.f <- ggplot(as.data.frame(seu_f@meta.data), aes(x=seurat_clusters, y=prediction.score))+geom_boxplot(outliers=F)+ylim(0,1)#+ggtitle("Maitra 2023: females")

p2.m <- ggplot(as.data.frame(seu_m@meta.data), aes(x=seurat_clusters, fill=seurat_low.res))+
	geom_bar(stat="count", position="stack")+
	scale_fill_manual(values=low.res.pal)+theme(legend.position="bottom", legend.key.size=unit(8,"pt"))
p2.f <- ggplot(as.data.frame(seu_f@meta.data), aes(x=seurat_clusters, fill=seurat_low.res))+
	geom_bar(stat="count", position="stack")+
	scale_fill_manual(values=low.res.pal)+theme(legend.position="bottom", legend.key.size=unit(8,"pt"))

# remove clusters with very low prediction score
cat("Remove low prediction score cluster(s), males:\n")
seu_m = seu_m[,seu_m$seurat_clusters!="16"]
dim(seu_m)
cat("Remove low	prediction score cluster(s), females:\n")
seu_f =	seu_f[,seu_f$seurat_clusters!="4"]
dim(seu_f)


p3.m <- ggplot(as.data.frame(seu_m@meta.data), aes(x=log10(nCount_RNA), group=Sample))+
	stat_ecdf()+facet_wrap(vars(seurat_low.res))+
	geom_hline(aes(yintercept=.5), color="red")+geom_vline(aes(xintercept=3), color="red")
p3.f <- ggplot(as.data.frame(seu_f@meta.data), aes(x=log10(nCount_RNA), group=Sample))+
        stat_ecdf()+facet_wrap(vars(seurat_low.res))+
        geom_hline(aes(yintercept=.5), color="red")+geom_vline(aes(xintercept=3), color="red")

ggsave(file="plots/08_compare_DEGs/Maitra/post-sn-qc_males.png", 
	gridExtra::grid.arrange(p0.m, p1.m, p2.m, p3.m, layout_matrix=matrix(c(1,2,3,4,4), ncol=1)),
	width=8, height=11)
ggsave(file="plots/08_compare_DEGs/Maitra/post-sn-qc_females.png", 
        gridExtra::grid.arrange(p0.f, p1.f, p2.f, p3.f, layout_matrix=matrix(c(1,2,3,4,4), ncol=1)),
        width=8, height=11)

#save(seu_m, file="processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT_post.Rdata")
#save(seu_f, file="processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT_post.Rdata")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
