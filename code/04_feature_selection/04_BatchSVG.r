setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(BatchSVG)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

#spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load any spe subet for rownames
#load("processed-data/04_feature_selection/per-sample_spe/V13B23-334_C1.Rdata")
#svg.set = rownames(tmp)
#rm(tmp)

cat("\nBatch effects examined: sample_id, slide, seq, condition, sex\n\n")
#list_batch_df <- featureSelect(input = spe, batch_effect = c("sample_id","slide","seq","condition","sex"), VGs = svg.set)

#saveRDS(list_batch_df, "processed-data/04_feature_selection/list_BatchSVG_dfs.rda")
#cat("\n\nOutput saved to: processed-data/04_feature_selection/list_BatchSVG_dfs.rda\n")

#load saved results
resList <- readRDS("processed-data/04_feature_selection/list_BatchSVG_dfs.rda")

#compile into single dframe
batchVars = names(resList)
comb.df = do.call(rbind, lapply(batchVars, function(x) {
  tmp = resList[[x]]
  colnames(tmp) = gsub(x, "batch", colnames(tmp))
  tmp$batch_variable = x
  return(tmp)
}))

#deviance
#sample_id at 5
comb.df$dev_outlier = comb.df$nSD_dev_batch>=5
comb.df$dev_outlier = ifelse(comb.df$batch_variable %in% c("seq","condition","slide","sex"), FALSE, comb.df$dev_outlier)
#condition at none
#seq needs to be higher
#filter(comb.df, batch_variable=="seq") %>% slice_max(n=4, nSD_dev_batch) #20 will do
comb.df$dev_outlier = ifelse(comb.df$batch_variable=="seq" & comb.df$nSD_dev_batch>20, TRUE, comb.df$dev_outlier)
#try slide at 20 too
comb.df$dev_outlier = ifelse(comb.df$batch_variable=="slide" & comb.df$nSD_dev_batch>20, TRUE, comb.df$dev_outlier)
#also sex at 20
comb.df$dev_outlier = ifelse(comb.df$batch_variable=="sex" & comb.df$nSD_dev_batch>20, TRUE, comb.df$dev_outlier)

#rank
comb.df$rank_outlier = comb.df$nSD_rank_batch>=5
comb.df$rank_outlier = ifelse(comb.df$batch_variable %in% c("seq","condition"), FALSE, comb.df$rank_outlier)
#seq needs to be higher
#filter(comb.df, batch_variable=="seq") %>% slice_max(n=4, nSD_rank_batch) #20
comb.df$rank_outlier = ifelse(comb.df$batch_variable=="seq" & comb.df$nSD_rank_batch>20, TRUE, comb.df$rank_outlier)
#same with condition
#filter(comb.df, batch_variable=="condition") %>% slice_max(n=4, nSD_rank_batch) #20
comb.df$rank_outlier = ifelse(comb.df$batch_variable=="condition" & comb.df$nSD_rank_batch>20, TRUE, comb.df$rank_outlier)


thold = as.data.frame(list("batch_label"=factor(c("sample_id","slide","seq","sex","condition"), 
                                                   levels=c("sample_id","slide","seq","sex","condition"),
                                                   labels=c("Sample", "Slide", "Sequencing core", "Sex", "Diagnosis")),
                           "dev_thold"=c(5,20,20,20,"none"),
                           "rank_thold"=c(5,5,20,5,20),
                           "dev_x"=rep(1e05,5),
                           "dev_y"=rep(5e06,5),
                           "rank_x"=rep(6000,5),
                           "rank_y"=rep(500,5)))

comb.df$batch_label = factor(comb.df$batch_variable, 
                                levels=c("sample_id","slide","seq","sex","condition"),
                                labels=c("Sample", "Slide", "Sequencing core", "Sex", "Diagnosis"))

p1 <- ggplot(comb.df, aes(x=dev_default, y=dev_batch, color=nSD_dev_batch>=5))+
  geom_point(pch=1)+scale_color_manual(values=c("grey","red3"))+
  geom_point(data=filter(comb.df, nSD_dev_batch>=5), pch=1)+
  geom_point(data=filter(comb.df, dev_outlier==TRUE), show.legend = FALSE)+
  geom_text(data=thold, aes(x=dev_x, y=dev_y, label=paste0("nSD cutoff= ",dev_thold,"\n(filled circle)")),
            hjust=0, vjust=1, color="black", size=3.5, lineheight=1)+
  scale_y_log10()+scale_x_log10()+
  labs(title="Deviance-based assessment", x="deviance (no block)", y="deviance (with block)",
       color="nSD diff. in dev. >5")+
  facet_wrap(vars(batch_label), ncol=1)+theme_bw()+
  theme(legend.position="bottom", strip.background=element_rect(fill="transparent", color="transparent"),
        panel.grid.minor=element_blank(), strip.text=element_text(size=12))


p2 <- ggplot(comb.df, aes(x=rank_default, y=rank_batch, color=nSD_rank_batch>=5))+
  geom_point(pch=1)+scale_color_manual(values=c("grey","red3"))+
  geom_point(data=filter(comb.df, nSD_rank_batch>=5), pch=1)+
  geom_point(data=filter(comb.df, rank_outlier==TRUE), show.legend = FALSE)+
  geom_text(data=thold, aes(x=rank_x, y=rank_y, label=paste0("nSD cutoff= ",rank_thold,"\n(filled circle)")),
            hjust=1, vjust=1, color="black", size=3.5, lineheight=1)+
  scale_y_reverse()+
  labs(title="Rank-based assessment", x="rank of deviance (no block)", y="rank of deviance (with block)",
       color="nSD diff. in rank >5")+
  facet_wrap(vars(batch_label), ncol=1)+theme_bw()+
  theme(legend.position="bottom", strip.background=element_rect(fill="transparent", color="transparent"),
        panel.grid.minor=element_blank(), strip.text=element_text(size=12))

ggsave(filename="plots/04_feature_selection/BatchSVG_nSD-cutoff.png", 
       gridExtra::grid.arrange(p1, p2, ncol=2), bg="white", height=12, width=8)
cat("\nnSD cutoff plot saved to: plots/04_feature_selection/BatchSVG_nSD-cutoff.png\n")

names(batchVars) = batchVars
batch.genes = do.call(c, lapply(batchVars, function(x) {
  tmp = filter(comb.df, batch_variable==x)
  out1 = tmp$gene_id[tmp$rank_outlier]
  names(out1) = tmp$gene_name[tmp$rank_outlier]
  out2 = tmp$gene_id[tmp$dev_outlier]
  names(out2) = tmp$gene_name[tmp$dev_outlier]
  list("rank"=out1, "dev"=out2)
})
)

saveRDS(batch.genes, "processed-data/04_feature_selection/BatchSVG_sample-slide-seq-sex-condition_list.rds")
cat("\nList of batch genes saved to: processed-data/04_feature_selection/BatchSVG_sample-slide-seq-sex-condition_list.rds\n")

tmp = filter(comb.df, gene_id %in% unlist(batch.genes))
default.df = distinct(tmp[,c("gene_name","gene_id","dev_default","rank_default")]) %>%
  arrange(rank_default)

p3 <- ggplot(tmp, aes(y=factor(gene_name, levels=rev(default.df$gene_name)), 
                x=dev_batch, color=batch_label))+
  geom_point(data=default.df, aes(x=dev_default), color="black", pch=4, size=2)+
  geom_point(pch=1, size=2)+scale_color_brewer(palette="Dark2")+
  scale_x_log10()+labs(x="deviance", color="block variable")+
  theme_bw()+theme(axis.text.y=element_text(face="italic"), axis.title.y=element_blank(),
                   legend.position="inside", legend.position.inside = c(.75,.1))

p4 <- ggplot(tmp, aes(y=factor(gene_name, levels=rev(default.df$gene_name)), 
                x=rank_batch, color=batch_label))+
  geom_point(data=default.df, aes(x=rank_default), color="black", pch=4, size=2)+
  geom_point(pch=1, size=2)+scale_color_brewer(palette="Dark2")+
  scale_x_reverse()+labs(x="rank of deviance", color="block variable")+
  theme_bw()+theme(axis.text.y=element_text(face="italic"), axis.title.y=element_blank(),
                   legend.position="inside", legend.position.inside = c(.75,.1))

ggsave(filename="plots/04_feature_selection/BatchSVG_batch-genes.png", 
       gridExtra::grid.arrange(p3, p4, ncol=2), bg="white", height=12, width=8)
cat("\nPlot of batch genes saved to: plots/04_feature_selection/BatchSVG_batch-genes.png\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
