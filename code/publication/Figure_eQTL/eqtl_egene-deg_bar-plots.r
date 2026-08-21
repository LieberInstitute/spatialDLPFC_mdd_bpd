setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})

egene.df = read.csv("processed-data/11_eQTL_coloc/seurat/tables/map_significant_pairs.csv.gz")
seurat_labels = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
names(seurat_labels) = seurat_labels
egene.df$cluster = factor(egene.df$context, levels=seurat_labels)


l1 = lapply(seurat_labels, function(x) filter(egene.df, cluster==x, pair_provenance=="cis_and_independent") %>% 
	distinct(variant_id, gene_name) %>% transmute(eGene= paste(variant_id, gene_name)) %>% pull(eGene))

l2 = lapply(seurat_labels, function(x) filter(egene.df, cluster==x, DEG==1, pair_provenance=="cis_and_independent") %>% 
	distinct(variant_id, gene_name) %>% transmute(eGene= paste(variant_id, gene_name)) %>% pull(eGene))

l3 = lapply(seurat_labels, function(x) filter(egene.df, cluster==x, DEG==1, pair_provenance=="cis_and_independent") %>% 
	distinct(variant_id, gene_name) %>% pull(gene_name) %>% unique())


df1 = data.frame("context"=factor(c(names(l1), names(l2), names(l3)), levels=seurat_labels, 
                            labels=c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")),
           "n_y"=c(sapply(l1, length), sapply(l2, length), sapply(l3, length)),
           "set"=factor(c(rep("SNP-gene pairs", length(l1)), rep("SNP-DEG pairs", length(l2)), rep("eGenes", length(l3))), 
                        levels=c("SNP-gene pairs","SNP-DEG pairs","eGenes")))


p1 <- ggplot(filter(df1, set=="SNP-gene pairs"), aes(x=context, y=n_y))+
  geom_bar(stat="identity", fill="grey50")+
  labs(y="SNP-gene pairs", title="all cis_and_indpendent (with ties)")+
  theme_bw()+theme(text=element_text(size=6), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

p2 <- ggplot(filter(df1, set!="SNP-gene pairs"), aes(x=context, y=n_y, fill=set))+
  geom_bar(stat="identity", position="dodge")+
  scale_y_continuous("SNP-gene pairs", limits=c(0,100), sec.axis= sec_axis(trans= ~ ., name="eGenes"))+
  scale_fill_manual(values=c("grey50","black"))+
  labs(title="DEG associations only")+
  theme_bw()+theme(text=element_text(size=6), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                   legend.key.size = unit(8,"pt"), legend.position="inside", legend.position.inside = c(.2,.8))


ggsave(file="plots/publication/Figure_eQTL/snp-gene-pairs_egenes_bar-plots.pdf",
       grid.arrange(p1, p2, layout_matrix=rbind(c(1,1,2,2,2))),
       width=4, height=2)



## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
