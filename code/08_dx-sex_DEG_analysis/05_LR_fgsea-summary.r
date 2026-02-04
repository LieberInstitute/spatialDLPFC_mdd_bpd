args = commandArgs(TRUE)
x= args[[1]]
print(x)

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(enrichR)
	library(pheatmap)
	library(igraph)
	library(gridExtra)
})

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/08_dx-sex_DEG_analysis/fgsea_functions.r")
source("code/08_dx-sex_DEG_analysis/jc-igraph_functions.r")

cpList = readRDS("plots/colorPalettes.rds")

#if(!dir.exists(paste0("plots/08_dx-sex_DEG_analysis/", gsub("\\.","-", gsub("_","-",x))))) dir.create(paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x))))

(res_file = "smoothed-k9-1663")
#(res_file = "seurat-pc30")

lr.results <- read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", res_file, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

lr.degs <- read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", res_file,
                          "_dx-sex_degs-F-test-t-test.csv"))
wm.sig = filter(lr.results, group==unlist(strsplit(x, "_"))[[1]], sex==unlist(strsplit(x, "_"))[[2]], adj.P.Val<.05, gene_id %in% lr.degs$gene_id,
	cluster=="WM")
cat("\nNumber of L-R WM DEGs (padj<.05) for", paste0(x,":"), nrow(wm.sig),"\n")
#nrow(adj.sig)

gmt_dbl = c("Reactome")
#gmt_dbl = c("Reactome","GO-BP","GO-CC")

for(gmt_db in gmt_dbl) {
	cat(paste0("\n\n",gmt_db))
	#load fgsea results
	fgsea.results = readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/", res_file, "_", gmt_db, "_WM_LR-fgsea-list.rda"))
	res1 = fgsea.results[[x]]
	cat("\nNumber of sig. GSEA terms for", paste0(gmt_db,":"),nrow(res1),"\n")
	#nrow(res1)

	ledge = tibble::deframe(res1[,c("pathway","leadingEdge")])
	l1 = lapply(ledge, intersect, y= wm.sig$gene_name)
	res1$LR_in_leadingEdge = sapply(l1, length)
	res1$LR_in_leadingEdge2 = sapply(l1, paste, collapse="/")

	ledge.df = select(res1, term=pathway, NES, LR_in_leadingEdge, LR_in_leadingEdge2) %>%
		mutate(source=rep(unlist(strsplit(res_file, "-"))[[1]], nrow(res1)))

	if(gmt_db=="Reactome") {
		p.names = sapply(strsplit(ledge.df$term, " R-HSA-"), function(x) x[[1]])
	} 
	if(gmt_db %in% c("GO-BP","GO-CC")) {
		p.names = sapply(strsplit(ledge.df$term, " \\(GO:"), function(x) x[[1]])
	}
	p.names = sapply(p.names, function(x) {
		c1 = unlist(strwrap(x, width=75))
		if(length(c1)>1) {
			return(paste(c1[[1]],"[...]"))
		} else {
			return(c1[[1]])
		}
	})
	ledge.df$name = p.names

	tmp1 = filter(ledge.df, LR_in_leadingEdge>0) 
	cat("\nNumber of sig. GSEA terms with WM L-R in leading edge for", paste0(gmt_db,":"), nrow(tmp1),"\n")
	#nrow(tmp1)

	tmp2 = tidyr::separate_rows(tmp1, LR_in_leadingEdge2, sep="/")
	colnames(tmp2)[4] = "gene_name"

	write.csv(tmp2, file=paste0("processed-data/08_dx-sex_DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_fgsea_WM_filtered-LR-DEG-leadingeEdge.csv"), row.names=F)
	cat("\nFiltered GSEA data frame (expanded by genes) saved to:", paste0("processed-data/08_dx-sex_DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_fgsea_WM_filtered-LR-DEG-leadingeEdge.csv"), "\n")

	#make plots
	## heatmap
	tmp3 = group_by(tmp2, term, name, gene_name) %>% 
		tidyr::pivot_wider(names_from="gene_name", values_from="NES", values_fill=0)
	m3 = as.matrix(tmp3[,5:ncol(tmp3)])
	rownames(m3) <- tmp3$name
	#color limits for NES
	color.limits = ceiling(max(abs(m3)))
	if(color.limits-max(abs(m3))>.5) color.limits = color.limits-.5
	seq.len = color.limits*4
	pal.len = seq.len-1

	col1 = colorRampPalette(rev(RColorBrewer::brewer.pal(n=5, "RdBu")))(pal.len)
	col1[median(seq(1,length(col1)))] = "white"
	phm3 = pheatmap(t(m3), color=col1,
                breaks= seq(-color.limits, color.limits, length.out=seq.len),
                treeheight_col = 12, treeheight_row = 12, silent=T, border_color = "grey",
                fontsize_row=7, fontsize_col=7, angle_col=90,
                main="Fill: NES (0 if not sig.)")

	#jaccard
	n.df = tmp1
	n.df$name = sapply(n.df$name, function(x) paste(strwrap(x, width=30), collapse="\n"))
	jc.df = jaccardFromList(ledge)
	e.df= filter(jc.df, jaccard>0)
	e.df = filter(e.df, reference %in% n.df$term, query %in% n.df$term)

	outList = generateIGRAPH(node_df= n.df, edge_df= e.df, size_var="LR_in_leadingEdge", 
		text_title=paste(gsub("_"," ",x), gmt_db, res_file, "WM"))

	#term to gene jaccard
	tmp2 = left_join(tmp2, lr.degs[,c("gene_name","med_prop.spots.detected_WM")])
	tmp2$node_size = ifelse(tmp2$med_prop.spots.detected_WM>.05, 6, 3)

	outlist2 = term2GeneIGRAPH(tmp2, source=FALSE, node_size="node_size", layout_style="kk", color_by="term",
		text_title=paste(gsub("_"," ",x), gmt_db, res_file, "WM"))

	#save plots
	pdf(file=paste0("plots/08_dx-sex_DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_",
		res_file, "_", gmt_db, "_WM_LR-fgsea-summary.pdf"), width=12, height=12)
	grid.arrange(phm3[[4]], top=paste(gsub("_"," ",x), gmt_db))
	plot(simplify(outList[["igraph"]]), layout=outList[["layout"]], 
		edge.width=E(outList[["igraph"]])$jc*5, 
		vertex.label.family="sans", vertex.frame.width=0,
		main=outList[["title_text"]],
		sub=outList[["sub_text"]]
	)
	legend("bottomleft", legend=c("Depleted","Enriched"),
		pch=16, pt.cex=1, cex=.7,
		col=c("#CFEBF7","#FFC0B5"))
	plot(outlist2$igraph, layout=outlist2$layout, vertex.label.family="sans", #vertex.size=6, 
		edges.curved=T,
		vertex.frame.width=0,
		main=outlist2[["title_text"]],
		sub="Smaller genes have median zero proportion <=5% of WM spots.")
	legend("bottomleft", legend=c("Depleted","Enriched","Gene"),
		pch=16, pt.cex=1, cex=.7,
		col=c("#CFEBF7","#FFC0B5","grey85"))
	dev.off()
	cat("\nPlots saved to:", paste0("plots/08_dx-sex_DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_",
                res_file, "_", gmt_db, "_WM_LR-fgsea-summary.pdf"),"\n")
}

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
