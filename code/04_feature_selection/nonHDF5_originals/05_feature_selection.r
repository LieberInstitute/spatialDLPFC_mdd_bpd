setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(UpSetR)
	library(here)
})

bindev.2k = read.csv(here("processed-data","04_preprocessing","subject-biased_genes-2000.csv"))
bindev.3k = read.csv(here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"))

genes.2k = unique(filter(bindev.2k, rank_brain<=2000)$gene)
genes.3k = unique(filter(bindev.3k, rank_brain<=3000)$gene)

outlier.2k = unique(filter(bindev.2k, all.slide.outlier>5)$gene)
outlier.3k = unique(filter(bindev.3k, all.slide.outlier>5)$gene)

outlier.any = union(outlier.2k, outlier.3k)

l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[
	unlist(lapply(l1, function(x) {
		if(dir.exists(here("processed-data","04_preprocessing",x))) return(FALSE)
		else {
			split_1 = unlist(strsplit(x, split="_"))[[1]]
			split.2 = unlist(strsplit(x, split="\\."))[[2]]
			if(split_1=="nnSVG" & split.2=="csv") return(TRUE)
			else {return(FALSE)}
		}
	}))
]

svg.df = do.call(rbind, lapply(l1, function(x) mutate(read.csv(here("processed-data","04_preprocessing",x)), file=x) %>% filter(padj<.05)))
svg.genes = unique(svg.df$gene_id)

lm = read.csv(here("processed-data","04_preprocessing","EXT_layer-markers_tableS5.csv"))
layer.markers = unique(filter(lm, rank<=3000)$ensembl)

name.key = union(distinct(bindev.2k, gene, gene_name), distinct(bindev.3k, gene, gene_name)) %>% 
	union(transmute(lm, gene=ensembl, gene_name=gene_name)) %>% 
	union(transmute(svg.df, gene=gene_id, gene_name=gene_name))

exclude.list <- list("MTRN"=name.key$gene[grep("^MTRN",name.key$gene_name)],
	"DNAJ"=name.key$gene[grep("^DNAJ",name.key$gene_name)],
	"HSP"=name.key$gene[grep("^HSP",name.key$gene_name)],
	"ribo"=name.key$gene[grep("RPS|RPL",name.key$gene_name)])


gene.list = list("bindev.2k"=genes.2k, "bindev.3k"=genes.3k, "svg"=svg.genes, "layer.markers"=layer.markers,"biased.any"=outlier.any,"exclude"=unique(unlist(exclude.list)))

#pdf(file=here("plots","04_preprocessing","feature-list_pre-filter_upset.pdf"), width=8, height=4)
#	upset(fromList(gene.list), nsets=6, nintersects=100, order.by="freq", mb.ratio=c(.6,.4), text.scale=2)
#dev.off()

gene.list2 = list("bindev.2k"=genes.2k, "bindev.3k"=genes.3k, "svg"=svg.genes)
#gene.list2 = lapply(gene.list2, setdiff, y=outlier.any)
gene.list2[['bindev.2k']]  = setdiff(gene.list2[['bindev.2k']], outlier.2k)
gene.list2[['bindev.3k']]  = setdiff(gene.list2[['bindev.3k']],	outlier.3k)
gene.list2 = lapply(gene.list2, setdiff, unlist(exclude.list))
gene.list2 = lapply(gene.list2, union, y=layer.markers)

gene.list2 = lapply(gene.list2, function(x) {
  tmp = filter(name.key, gene %in% x)
  x1 = tmp$gene
  names(x1) = tmp$gene_name
  return(x1)
})

gene.list3 = c(gene.list2, list("layer.markers"=layer.markers, "biased.any"=outlier.any, "exclude"=unique(unlist(exclude.list))))
pdf(file=here("plots","04_preprocessing","feature-list_pre-and-post-filter_upset.pdf"), width=8, height=4)
        upset(fromList(gene.list), nsets=6, nintersects=100, order.by="freq", mb.ratio=c(.6,.4), text.scale=2)
	upset(fromList(gene.list3), nsets=6, nintersects=100, order.by="freq", mb.ratio=c(.6,.4), text.scale=2)
dev.off()

cat("\n\nFinal feature list size:\n")
unlist(lapply(gene.list2, length))
saveRDS(gene.list2, here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))

gene.list4 = list("biased.2k"=outlier.2k, "biased.3k"=outlier.3k, "exclude"=exclude.list)
saveRDS(gene.list4, here("processed-data","04_preprocessing","biased_excluded_feature-list.rda"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
