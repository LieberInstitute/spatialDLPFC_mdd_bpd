setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(here)
})

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
bindev.3k = read.csv(here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"), row.names=1)

plot.genes = unique(filter(bindev.3k, all.slide.outlier>5)$gene)
names(plot.genes) = rowData(spe)[plot.genes,"gene_name"]

best.rank.slide.df = filter(bindev.3k, gene %in% plot.genes) %>% group_by(slide, gene, gene_name) %>% summarize(best.rank.slide=min(rank_brain), .groups="drop")

best.rank.all.df = filter(bindev.3k, gene %in% plot.genes) %>% group_by(gene, gene_name) %>% summarize(best.rank.all = min(rank_brain), .groups="drop") %>% 
	arrange(best.rank.all) %>% mutate(index=row_number(), i2=length(plot.genes)-index, ytext=paste(best.rank.all,gene_name, sep=" - ")) %>% 
	arrange(i2) %>% mutate(ylabel=factor(i2, levels=i2,labels=ytext))

plot.rank.df = filter(bindev.3k, gene %in% plot.genes) %>% select(slide, gene, gene_name, rank_brain, rank_default, r.diff) %>% 
	left_join(best.rank.all.df, by=c("gene","gene_name"))

p1, ggplot(plot.rank.df, aes(x=slide, y=ylabel, size=rank_default, color=rank_brain))+
	geom_count()+scale_color_viridis_c(direction=-1)+
	scale_size(breaks=c(100,500,1000,1500,2000))+
	labs(x="slide",y="(best rank across slides) - gene name", size="rank (no batch)", color="rank (batch = sample)",
		subtitle="missing spot indicates gene not in top 3k deviant genes for default or subject-batch")+
	theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

l2 = unique(spe$sample_id)
names(l2) = lapply(l2, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

plot.genes.df = do.call(rbind, lapply(seq_along(l2), function(x) {
	m1 = as.data.frame(do.call(rbind, 
	lapply(plot.genes, function(y) {
		gex = logcounts(l2[[x]])[y,]
		c("avg"=mean(gex), "n"=sum(gex>0)/length(gex))
		})
	))
	m1$gene_name = names(plot.genes)
	m1$gene = plot.genes
	m1$brain = names(l2)[x]
	m1$slide = unique(l2[[x]]$slide)
	m1$position = unique(l2[[x]]$position)
	m1$condition = unique(l2[[x]]$condition)
	return(m1)
}))

plot.genes.df <- left_join(plot.genes.df, best.rank.slide.df, by=c("slide","gene","gene_name")) %>% 
	left_join(best.rank.all.df, by=c("gene","gene_name")) %>% 
	mutate(xlabel=paste(position, brain, condition))

p2 <- ggplot(plot.genes.df, aes(x=xlabel, y=ylabel, size=n, color=avg))+
	geom_count()+scale_color_viridis_c(option="F", direction=-1)+
	facet_wrap(vars(slide), ncol=6, scales="free_x")+
	labs(x="", y="(best rank across slides) - gene name",size="expr (prop. spots)", color="expr (avg. logcounts)")+
	theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file=here("plots","04_preprocessing","biased-features_dotplots.png"), plot=gridExtra::grid.arrange(p1,p2,ncol=2), width=20, height=12, bg="white")

cat("plot destination:", here("plots","04_preprocessing","biased-features_dotplots.png"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
