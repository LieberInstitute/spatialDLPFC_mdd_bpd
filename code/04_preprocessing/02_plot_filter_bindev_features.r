setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(SpatialExperiment)
	library(here)
})

source(here("code","04_preprocessing","02_plot_bindev_functions.r"))
l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[grep("^bindev_V.{9}_default-brain",l1)]

# perform for top 3k features
cat("Find subject-biased genes in top 3k deviant features\n")
bindev.3k = do.call(rbind, lapply(l1, function(x) {
	tmp = read.csv(here("processed-data","04_preprocessing",x))
	y=substr(x,8,17)
	tmp = filter(tmp, rank_default<=3000 | rank_brain<=3000)
	tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
	return(tmp)
})
)

f1 = plotBiasedFeatures(bindev.3k, 3000)
ggsave(file=here("plots","04_preprocessing","subject-biased_genes_hist-scatter-3000.png"), plot=f1[[1]], height=10, width=12, bg="white")

write.csv(f1[[2]], here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"))

# perform for top 2k features
cat("\nFind subject-biased genes in top 2k deviant features\n")
bindev.2k = do.call(rbind, lapply(l1, function(x) {
        tmp = read.csv(here("processed-data","04_preprocessing",x))
        y=substr(x,8,17)
        tmp = filter(tmp, rank_default<=2000 | rank_brain<=2000)
        tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
        return(tmp)
})
)

f2 = plotBiasedFeatures(bindev.2k, 2000)
ggsave(file=here("plots","04_preprocessing","subject-biased_genes_hist-scatter-2000.png"), plot=f2[[1]], height=10, width=12, bg="white")

write.csv(f2[[2]], here("processed-data","04_preprocessing","subject-biased_genes-2000.csv"))

### below is code for dotplot that i am temporarily omitting
#load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
#l2 = unique(spe$sample_id)
#names(l2) = lapply(l2, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
#l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

#plot.genes = filter(d1, rank_brain<=3000) %>% pull(gene) %>% unique()
#names(plot.genes) = rowData(spe)[plot.genes,"gene_name"]

#plot.genes.df = do.call(rbind, lapply(seq_along(l2), function(x) {
#	m1 = as.data.frame(do.call(rbind, 
#		lapply(plot.genes, function(y) {
#			gex = logcounts(l2[[x]])[y,]
#			c("avg"=mean(gex), "n"=sum(gex>0)/length(gex))
#		})
#	))
#	m1$gene_name = names(plot.genes)
#	m1$gene = plot.genes
#	m1$brain = names(l2)[x]
#	m1$slide = unique(l2[[x]]$slide)
#	return(m1)
#}))

#plot1 <- ggplot(plot.genes.df, aes(x=brain, y=gene_name, size=n, color=avg))+
#	geom_count()+scale_color_viridis_c(option="F", direction=-1)+
#	facet_grid(cols=vars(slide), scales="free_x")+
#	labs(title="Biased genes with top 3k rank in any slide")+
#	theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

#plot.genes2 = filter(d1, rank_brain>=3000) %>% pull(gene) %>% unique()
#plot.genes = setdiff(plot.genes2, plot.genes)
#names(plot.genes) = rowData(spe)[plot.genes,"gene_name"]

#plot.genes.df = do.call(rbind, lapply(seq_along(l2), function(x) {
#        m1 = as.data.frame(do.call(rbind, 
#                lapply(plot.genes, function(y) {
#                        gex = logcounts(l2[[x]])[y,]
#                        c("avg"=mean(gex), "n"=sum(gex>0)/length(gex))
#                })
#        ))
#        m1$gene_name = names(plot.genes)
#        m1$gene = plot.genes
#        m1$brain = names(l2)[x]
#        m1$slide = unique(l2[[x]]$slide)
#        return(m1)
#}))

#plot2 <- ggplot(plot.genes.df, aes(x=brain, y=gene_name, size=n, color=avg))+
#        geom_count()+scale_color_viridis_c(option="F", direction=-1)+
#        facet_grid(cols=vars(slide), scales="free_x")+
#        labs(title="Biased genes not top 3k  rank in any slide")+
#        theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

#pdf(here("plots","04_preprocessing","subject-biased_genes_dotplot-3000.pdf"), height=8, width=7)
#	plot1
#	plot2
#dev.off()

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
