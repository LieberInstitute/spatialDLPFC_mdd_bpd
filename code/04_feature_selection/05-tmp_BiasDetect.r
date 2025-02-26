setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(BiasDetect)
	library(UpSetR)
})
set.seed(123)

dfList = list(dummy.slide=read.csv("processed-data/04_feature_selection/test_bindev_batch-dummy-slide.csv"),
              sample=read.csv("processed-data/04_feature_selection/test_bindev_batch-sample.csv"),
	      seq=read.csv("processed-data/04_feature_selection/test_bindev_batch-seq-rnd.csv"),
              sex=read.csv("processed-data/04_feature_selection/test_bindev_batch-sex.csv"),
              condition=read.csv("processed-data/04_feature_selection/test_bindev_batch-condition.csv"))
dfList = lapply(dfList, function(batch_df) {
  batch_df$d_diff <- (batch_df$dev_default-batch_df$dev_batch)/
    batch_df$dev_batch
  mean_dev <- mean(batch_df$d_diff)
  sd_dev <- sd(batch_df$d_diff)
  batch_df$nSD_dev <- (batch_df$d_diff - mean_dev) / sd_dev
  
  batch_df$r_diff <- batch_df$rank_batch-batch_df$rank_default
  mean_rank <- mean(batch_df$r_diff)
  sd_rank <- sd(batch_df$r_diff)
  batch_df$nSD_rank <- (batch_df$r_diff - mean_rank) / sd_rank
  
  colnames(batch_df)[1] = "gene"
  
  return(batch_df)
})

#plots
plotList = lapply(dfList, biasDetect, nSD_dev=5, visual=T)
ggsave("plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_scatter.png",
       gridExtra::grid.arrange(plotList$dummy.slide[[1]]+ggtitle("dummy.slide: deviance"), plotList$sample[[1]]+ggtitle("sample: deviance"), 
			       plotList$seq[[1]]+ggtitle("seq rnd: deviance"),
                               plotList$sex[[1]]+ggtitle("sex: deviance"), plotList$condition[[1]]+ggtitle("condition: deviance"),
                               plotList$dummy.slide[[2]]+ggtitle("dummy.slide: rank"), plotList$sample[[2]]+ggtitle("sample: rank"),
			       plotList$seq[[2]]+ggtitle("seq rnd: rank"),
                               plotList$sex[[2]]+ggtitle("sex: rank"), plotList$condition[[2]]+ggtitle("condition: rank"),
                               ncol=5),
       bg="white", height=6, width=20)
cat("\n\nScatter plot with nSD bins saved to: plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_scatter.png\n")

#nSD filters based on scatter plots
geneList = list()
cat("\n\nDummy Slide: nSD_dev>=5 and nSD_rank>=5\n")
geneList$dummy.slide_dev = filter(dfList[[1]], nSD_dev>=5)$gene
names(geneList$dummy.slide_dev) <- filter(dfList[[1]], nSD_dev>=5)$gene_name
geneList$dummy.slide_rank = filter(dfList[[1]], nSD_rank>=5)$gene
names(geneList$dummy.slide_rank) = filter(dfList[[1]], nSD_rank>=5)$gene_name

cat("\nSample: nSD_dev>=5 and nSD_rank>=5\n")
geneList$sample_dev = filter(dfList[[2]], nSD_dev>=5)$gene
names(geneList$sample_dev) <- filter(dfList[[2]], nSD_dev>=5)$gene_name
geneList$sample_rank = filter(dfList[[2]], nSD_rank>=5)$gene
names(geneList$sample_rank) = filter(dfList[[2]], nSD_rank>=5)$gene_name

cat("\nSeq. round: nSD_dev>=6 and nSD_rank>=6\n")
geneList$seq_dev = filter(dfList[[3]], nSD_dev>=6)$gene
names(geneList$seq_dev) <- filter(dfList[[3]], nSD_dev>=6)$gene_name
geneList$seq_rank = filter(dfList[[3]], nSD_rank>=6)$gene
names(geneList$seq_rank) = filter(dfList[[3]], nSD_rank>=6)$gene_name

cat("\nSex: nSD_dev>=10 and nSD_rank>=5\n")
geneList$sex_dev = filter(dfList[[4]], nSD_dev>=10)$gene
names(geneList$sex_dev) <- filter(dfList[[4]], nSD_dev>=10)$gene_name
geneList$sex_rank = filter(dfList[[4]], nSD_rank>=5)$gene
names(geneList$sex_rank) = filter(dfList[[4]], nSD_rank>=5)$gene_name

#condition hardly looks like it deserves any filters. check to see what genes are unique
cat("\nCondition: nSD_dev>=25 and nSD_rank>=15\n")
cat("*** Condition doesn't seem to have any unique effect")
geneList$condition_dev = filter(dfList[[5]], nSD_dev>=25)$gene
names(geneList$condition_dev) <- filter(dfList[[5]], nSD_dev>=25)$gene_name
geneList$condition_rank = filter(dfList[[5]], nSD_rank>=15)$gene
names(geneList$condition_rank) = filter(dfList[[5]], nSD_rank>=15)$gene_name

cat("\n\nLength of batch effect genes per block:\n")
sapply(geneList, length)

saveRDS(geneList, "processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-seq-sex-condition_list.rds")
cat("\nGene list saved to: processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-seq-sex-condition_list.rds\n")

png(filename="plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_nSD-filtered_upset.png",
	bg="white", width=7, height=6, unit="in", res=300)
upset(fromList(geneList), text.scale = 2, mb.ratio=c(.6,.4), sets=names(geneList)[10:1], keep.order = T)
dev.off()
cat("\nUpset plot saved to: plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_nSD-filtered_upset.png\n")

#postfilter plots
p1 = ggplot(dfList[["dummy.slide"]], aes(x=dev_default, y=dev_batch, color=gene %in% geneList$dummy.slide_dev))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_x_log10() + scale_y_log10()+
  #ggrepel::geom_text_repel(data=filter(dfList[["dummy.slide"]], gene %in% geneList$dummy.slide_dev),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = 1, intercept = 0), lty = 2)+
  labs(x= "dev (no batch)", y="dev (batch)", title="dummy.slide: deviance (nSD>=5)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))
p2 = ggplot(dfList[["dummy.slide"]], aes(x=rank_default, y=rank_batch, color=gene %in% geneList$dummy.slide_rank))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_y_reverse()+
  #ggrepel::geom_text_repel(data=filter(dfList[["dummy.slide"]], gene %in% geneList$dummy.slide_rank),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = -1, intercept = 0), lty = 2)+
  labs(x= "rank (no batch)", y="rank (batch)", title="dummy.slide: rank (nSD>=5)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))

p3 = ggplot(dfList[["sample"]], aes(x=dev_default, y=dev_batch, color=gene %in% geneList$sample_dev))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_x_log10() + scale_y_log10()+
  #ggrepel::geom_text_repel(data=filter(dfList[["sample"]], gene %in% geneList$sample_dev),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = 1, intercept = 0), lty = 2)+
  labs(x= "dev (no batch)", y="dev (batch)", title="sample: deviance (nSD>=5)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))
p4 = ggplot(dfList[["sample"]], aes(x=rank_default, y=rank_batch, color=gene %in% geneList$sample_rank))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_y_reverse()+
  #ggrepel::geom_text_repel(data=filter(dfList[["sample"]], gene %in% geneList$sample_rank),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = -1, intercept = 0), lty = 2)+
  labs(x= "rank (no batch)", y="rank (batch)", title="sample: rank (nSD>=5)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))


p5 = ggplot(dfList[["seq"]], aes(x=dev_default, y=dev_batch, color=gene %in% geneList$seq_dev))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_x_log10() + scale_y_log10()+
  #ggrepel::geom_text_repel(data=filter(dfList[["seq"]], gene %in% geneList$seq_dev),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = 1, intercept = 0), lty = 2)+
  labs(x= "dev (no batch)", y="dev (batch)", title="seq: deviance (nSD>=6)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))
p6 = ggplot(dfList[["seq"]], aes(x=rank_default, y=rank_batch, color=gene %in% geneList$seq_rank))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_y_reverse()+
  #ggrepel::geom_text_repel(data=filter(dfList[["seq"]], gene %in% geneList$seq_rank),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = -1, intercept = 0), lty = 2)+
  labs(x= "rank (no batch)", y="rank (batch)", title="seq: rank (nSD>=6)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))


p7 = ggplot(dfList[["sex"]], aes(x=dev_default, y=dev_batch, color=gene %in% geneList$sex_dev))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_x_log10() + scale_y_log10()+
  #ggrepel::geom_text_repel(data=filter(dfList[["sex"]], gene %in% geneList$sex_dev),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = 1, intercept = 0), lty = 2)+
  labs(x= "dev (no batch)", y="dev (batch)", title="sex: deviance (nSD>=10)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))
p8 = ggplot(dfList[["sex"]], aes(x=rank_default, y=rank_batch, color=gene %in% geneList$sex_rank))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_y_reverse()+
  #ggrepel::geom_text_repel(data=filter(dfList[["sex"]], gene %in% geneList$sex_rank),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = -1, intercept = 0), lty = 2)+
  labs(x= "rank (no batch)", y="rank (batch)", title="sex: rank (nSD>=5)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))

p9 = ggplot(dfList[["condition"]], aes(x=dev_default, y=dev_batch, color=gene %in% geneList$condition_dev))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_x_log10() + scale_y_log10()+
  #ggrepel::geom_text_repel(data=filter(dfList[["condition"]], gene %in% geneList$condition_dev),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = 1, intercept = 0), lty = 2)+
  labs(x= "dev (no batch)", y="dev (batch)", title="condition: deviance (nSD>=25)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))
p10 = ggplot(dfList[["condition"]], aes(x=rank_default, y=rank_batch, color=gene %in% geneList$condition_rank))+
  geom_point(size=.5)+scale_color_manual(values=c("grey50","red3"))+
  scale_y_reverse()+
  #ggrepel::geom_text_repel(data=filter(dfList[["condition"]], gene %in% geneList$condition_rank),
  #                         aes(label = gene_name), size = 3)+
  geom_abline(aes(slope = -1, intercept = 0), lty = 2)+
  labs(x= "rank (no batch)", y="rank (batch)", title="condition: rank (nSD>=15)")+
  theme_bw()+theme(legend.position="none", plot.title=element_text(size=12))

ggsave("plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_nSD-filtered_scatter.png",
	gridExtra::grid.arrange(p1, p3, p5, p7, p9, p2, p4, p6, p8, p10, ncol=5),
	bg="white", height=6, width=16)
cat("\nScatter plot of nSD filtered batch effect genes saved to: plots/04_feature_selection/bindev-batch-effect_dummyslide-sample-seq-sex-condition_nSD-filtered_scatter.png\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
