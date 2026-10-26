setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

set.seed(123)

### load in F-test sig genes to test subset
#ftest = read.table("raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt")

## or just domain-CT F-test
ftest = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test_all.csv") %>% filter(adj.P.Val<.05)
colnames(ftest)[56] = "V1"

mbv = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_moderated-t-test.csv")
mbv$cluster = factor(mbv$cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
	labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))
mbv$source = "MBv"

mt23 = read.csv("processed-data/08_compare_DEGs/Maitra/layer-restricted-pc3-ncells_low-res_moderated-t-test.csv")
mt23$cluster_og = factor(mt23$cluster, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"))
mt23$cluster_match = factor(mt23$cluster, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"),
	c("M.V","Ast","L2.3","L2.3","L4","Inb","L5","L6","Olg"))
mt23$source = "Maitra 2023"

rz24 = read.csv("processed-data/08_compare_DEGs/SZBDMulti-seq/layer-restricted-pc9-age-ncells_low-res_moderated-t-test.csv")
rz24$cluster_og = factor(rz24$cluster, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"))
rz24$cluster_match = factor(rz24$cluster, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"),
	c("M.V","Ast","L2.3","L2.3","L4","Inb","L5","L6","Olg"))
rz24$source = "SZBDMulti-seq"


mdd = inner_join(mbv[,c("gene_id","gene_name","cluster","sex","group","logFC","t","P.Value")],
	mt23[,c("gene_id","cluster_og","cluster_match","sex","group","logFC","t","P.Value","source")],
	by=c("gene_id","cluster"="cluster_match","sex","group"), suffix=c("_mbv","_sn"))

tmp = filter(mdd, gene_name %in% ftest$V1)
tmp1 = filter(tmp, sex=="F")
tmp2 = filter(tmp, sex=="M")
cor.df1 = do.call(rbind, lapply(levels(tmp$cluster_og), function(x) {
	tt1 = filter(tmp1, cluster_og==x)
	tt2 = filter(tmp2, cluster_og==x)
	as.data.frame(list("source"=c("Maitra 2023","Maitra 2023"),
		"sex"=c("F","M"),
		"cluster_og"=c(x,x),
		"corr_t"=c(cor(tt1$t_mbv, tt1$t_sn),
			cor(tt2$t_mbv, tt2$t_sn))
	))
}))


bd = inner_join(mbv[,c("gene_id","gene_name","cluster","sex","group","logFC","t","P.Value")],
	rz24[,c("gene_id","cluster_og","cluster_match","sex","group","logFC","t","P.Value","source")],
	by=c("gene_id","cluster"="cluster_match","sex","group"), suffix=c("_mbv","_sn"))


tmp = filter(bd, gene_name %in% ftest$V1)
tmp1 = filter(tmp, sex=="F")
tmp2 = filter(tmp, sex=="M")
cor.df2 = do.call(rbind, lapply(levels(tmp$cluster_og), function(x) {
	tt1 = filter(tmp1, cluster_og==x)
	tt2 = filter(tmp2, cluster_og==x)
	as.data.frame(list("source"=c("SZBDMulti-seq","SZBDMulti-seq"),
		"sex"=c("F","M"),
		"cluster_og"=c(x,x),
		"corr_t"=c(cor(tt1$t_mbv, tt1$t_sn),
			cor(tt2$t_mbv, tt2$t_sn))
	))
}))


all.df = bind_rows(mdd, bd) %>% filter(gene_name %in% ftest$V1) %>%
	mutate(cluster_og = factor(cluster_og, levels=levels(mt23$cluster_og)))
cor.df = bind_rows(cor.df1, cor.df2) %>%
	mutate(cluster_og = factor(cluster_og, levels=levels(mt23$cluster_og)))

p3 <- ggplot(all.df, aes(x=t_mbv, y=t_sn))+
        geom_point(size=.3)+
        geom_smooth(method="lm", se=F, color="red3")+
        geom_text(data=cor.df, aes(x=-7, y=7, label=paste0("rho: ", round(corr_t, 3))), hjust=0, size=3, color="red3")+
        facet_grid(cols=vars(paste(source, sex)), rows=vars(cluster_og))+
	ylim(-8,8)+xlim(-8,8)+
        ggtitle("Domain-restircted model DE compared to sex-specific snRNA-seq")+
        theme_bw()+theme(panel.grid.minor=element_blank(), text=element_text(size=10))


ggsave(file="plots/08_compare_DEGs/snRNAseq-sex-stratified_domain-restricted-model_domain-ct-F-sig-genes.png",
	p3, height= 8, width= 5.5, bg="white")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
