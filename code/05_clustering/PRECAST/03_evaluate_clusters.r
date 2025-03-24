setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
	library(ggspavis)
	library(dplyr)
})

#loss plot
#x = readLines("code/05_clustering/PRECAST/logs/precast_n1079_k9_15278541.log")
#x = readLines("code/05_clustering/PRECAST/logs/precast_n1663_k9_15309880.log")
x = readLines("code/05_clustering/PRECAST/logs/precast_n1629-qual-no-batch_k9_15326871.log")
iters = x[grep("iter =", x)]
iter.list = strsplit(iters, ", ")
df = data.frame(iter=2:(length(iters)+1), 
           loglik=as.numeric(sapply(iter.list, function(x) substr(x[[2]], start=9, stop=20))),
           d.loglik=as.numeric(sapply(iter.list, function(x) substr(x[[3]], start=9, stop=16))))

#png(file="plots/05_clustering/PRECAST_n1079-k9_lglk-loss-plot.png", bg="white")
#	plot(df[["iter"]], df[["loglik"]], main="PRECAST n=1079 k=9", xlab="iteration", ylab="loglik")
#dev.off()
#png(file="plots/05_clustering/PRECAST_n1663-k9_lglk-loss-plot.png", bg="white")
#        plot(df[["iter"]], df[["loglik"]], main="PRECAST n=1663 k=9", xlab="iteration", ylab="loglik")
#dev.off()
png(file="plots/05_clustering/PRECAST_n1629-k9_lglk-loss-plot.png", bg="white")
        plot(df[["iter"]], df[["loglik"]], main="PRECAST n=1629 k=9", xlab="iteration", ylab="loglik")
dev.off()
#load seurat results
#load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1079.Rdata")
#load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1629.Rdata")

#load spe and add seurat key
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spe$seurat_key = paste(spe$sample_id, colnames(spe), sep="_")
#match obs
mdata = seuInt@meta.data[spe$seurat_key,]
stopifnot(identical(rownames(mdata), spe$seurat_key))
#spe$precast_k9_1079 = as.factor(mdata$cluster)
#spe$precast_k9_1663 = as.factor(mdata$cluster)
spe$precast_k9_1629 = as.factor(mdata$cluster)

quickResaveHDF5SummarizedExperiment(spe)
#cat("\n\nPRECAST clusters (n=1079 genes, k=9) updated to spe with quickResave\n\n")
#cat("\n\nPRECAST clusters (n=1663 genes, k=9) updated to spe with quickResave\n\n")
cat("\n\nPRECAST clusters (n=1629 genes, k=9) updated to spe with quickResave\n\n")

#table(colData(spe)[,c("precast_k9_1079","problem_area_flag")])
#table(colData(spe)[,c("precast_k9_1663","problem_area_flag")])
table(colData(spe)[,c("precast_k9_1629","problem_area_flag")])

#plot it out bb!!
#uniquepal = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#8B0000","#FF7F00","#CAB2D6")
#names(uniquepal) = c("2","7","9","5","1","4","6","3","8")
#uniquepal = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FF7F00","#FDBF6F","#CAB2D6")
#names(uniquepal) = c("7","3","1","2","6","9","5","4","8")
uniquepal = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C","#FDBF6F")
names(uniquepal) = c("3","5","2","6","7","1","4","9","8")

spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "A1", spe$array)

seed = levels(as.factor(spe$dummy_slide))
slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList = lapply(slideList, function(x) {
	do.call(cbind, lapply(x, function(y) spe[,spe$dummy_slide==y]))
})

cat("\nGenerate spot plots...\n")
clusPlot = lapply(slideList, function(x) {
	suppressMessages(
		plotSpots(x, annotate="precast_k9_1629", 
			#annotate="precast_k9_1663",
			#annotate="precast_k9_1079", 
			point_size=.5, sample_id="sample_id")+
		scale_color_manual(values=uniquepal)+
		facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
		theme(panel.background=element_rect(fill="grey30"))
	)
})

pdf(file="plots/05_clustering/PRECAST_n1629-k9_spot-plots.pdf", width=12, height=16)
#pdf(file="plots/05_clustering/PRECAST_n1663-k9_spot-plots.pdf", width=12, height=16)
#pdf(file="plots/05_clustering/PRECAST_n1079-k9_spot-plots.pdf", width=12, height=16)
	clusPlot[[1]]
	clusPlot[[2]]
	clusPlot[[3]]
	clusPlot[[4]]
	clusPlot[[5]]
dev.off()
#cat("\nSaved spot plots to: plots/05_clustering/PRECAST_n1079-k9_spot-plots.pdf\n")
#cat("\nSaved spot plots to: plots/05_clustering/PRECAST_n1663-k9_spot-plots.pdf\n")
cat("\nSaved spot plots to: plots/05_clustering/PRECAST_n1629-k9_spot-plots.pdf\n")

#show work for initial layer annotation
#estimate spatial domain
layer.markers = read.csv("processed-data/04_feature_selection/EXT_TableS9_sig_genes_FDR5perc_enrichment.csv") %>%
  filter(stat>0, spatial_domain_resolution=="Sp09") %>%
  mutate(domain_simple=factor(test, 
                              levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
                              labels=c("L1 (1)","L1 (2)","L2","L3","L4","L5","L6","WM (1)","WM (2)")))

domains = c("L1 (1)","L1 (2)","L2","L3","L4","L5","L6","WM (1)","WM (2)")
names(domains) = domains
top100 = lapply(domains, function(x) filter(layer.markers, domain_simple==x) %>% slice_max(n=100, stat) %>% pull(gene))

#p1 <- UpSetR::upset(UpSetR::fromList(top100), nsets=9)
##which of the L1 clusters is more vascular
#grepl("VIM", top100) #L1 (1)
#pull out only unique markers
t1 = table(unlist(top100))
top100.unique = names(t1)[t1==1]
#length(top100.unique)
#add in WM genes that were in both WM groups
t2 = table(unlist(top100[c("WM (1)","WM (2)")]))
top100.wm = names(t2)[t2==2]
top100.unique = c(top100.unique, top100.wm)
#add GABA genes
gaba = c("SLC6A1","GAD1","GAD2","SST","CCK","TAC1","PVALB","CNR1","RELN")
#setdiff(gaba, rowData(spe)$gene_name)
#intersect(gaba, rowData(spe)$gene_name)
top100.unique = union(top100.unique, gaba)
cat("\n\nNumber of (mostly) unique layer markers for initial annotation:\n")
length(top100.unique)
#make sure that all genes are present in my dataset
#length(intersect(top100.unique, rowData(spe)$gene_name))==length(top100.unique)
#make df
top100.unique.df = filter(layer.markers, gene %in% top100.unique) %>% group_by(gene) %>% slice_max(n=1, stat) %>% select(gene, domain_simple, fdr, stat) %>%
  mutate(domain_simple=ifelse(gene %in% gaba, "GABA", as.character(domain_simple)))
#add gene_id
lut = rowData(spe)[rowData(spe)$gene_name %in% top100.unique,]
rownames(lut) = lut$gene_name
top100.unique.df$gene_id = lut[top100.unique.df$gene,"gene_id"]

#spe$precast_cluster = paste0("c",spe$precast_k9_1079)
#spe$precast_cluster = paste0("c",spe$precast_k9_1663)
spe$precast_cluster = paste0("c",spe$precast_k9_1629)
table(spe$precast_cluster)
cat("\n\nSummarise experiment to get mean expression for heatmap:\n")
Sys.time()
spe_summ = scuttle::aggregateAcrossCells(spe[top100.unique.df$gene_id,], ids=colData(spe)[,c("sample_id","precast_cluster")], 
                                         statistics=c("mean"),
                                         use.assay.type="logcounts")
#spe_summ
Sys.time()
prc_clus = paste0("c",1:9)
names(prc_clus) = prc_clus
m1 = sapply(prc_clus, function(x) {
  idx = colData(spe_summ)$precast_cluster==x
  rowMeans(logcounts(spe_summ)[,idx])
})

row_annot = data.frame("cluster"=top100.unique.df$domain_simple)
rownames(row_annot) = top100.unique.df$gene_id
annot_colors=list("cluster"=c("L1 (1)"="#FF7F00","L1 (2)"= "#1F78B4", "L2"="#33A02C", "L3"="#A6CEE3", "L4"="grey", "L5"="#B2DF8A", "L6"="#FB9A99", "WM (1)"="#8B0000", "WM (2)"="#E31A1C", "GABA"="#CAB2D6"))

p2 <- pheatmap::pheatmap(m1[rownames(row_annot)[order(row_annot[,1])],], annotation_row=row_annot, annotation_colors = annot_colors, 
                   show_rownames = F, cluster_rows=F, annotation_names_row = F,
                   scale="row", angle_col = 0, treeheight_col = 10)

#ggsave("plots/05_clustering/PRECAST_n1663-k9_initial-layer-annotation.pdf", p2[[4]], width=6, height=7, units="in")
#cat("\nInitial layer annotation heatmap saved to: plots/05_clustering/PRECAST_n1663-k9_initial-layer-annotation.pdf")
ggsave("plots/05_clustering/PRECAST_n1629-k9_initial-layer-annotation.pdf", p2[[4]], width=6, height=7, units="in")
cat("\nInitial layer annotation heatmap saved to: plots/05_clustering/PRECAST_n1629-k9_initial-layer-annotation.pdf")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
