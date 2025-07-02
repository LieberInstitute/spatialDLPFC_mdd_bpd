setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(ggplot2)
	library(dplyr)
	library(escheR)
	library(gridExtra)
	library(scater)
	library(bluster)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

cpList = readRDS("plots/colorPalettes.rds")
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(colnames(spe), rownames(cdata)))
spe$precast_k9_1663_f = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","L5","L6","WM","GABA","low UMI"))
table(spe$precast_k9_1663_f)

#check GAD1 and GAD2
spe$is_gaba = spe$precast_k9_1663_f=="GABA"
spe$gad1 = logcounts(spe)[rowData(spe)$gene_name=="GAD1",]
spe$gad2 = logcounts(spe)[rowData(spe)$gene_name=="GAD2",]

#gaba spot plots
sub_samples = paste0("V13B23-", c("329_A1","332_A1","309_D1","352_B1","352_D1","308_A1"))

plist <- do.call(c, lapply(sub_samples, function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="precast_k9_1663_f", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="is_gaba", point_size = .5)
  p1 <- p+scale_color_manual("", values=c(cpList$earthy.pal1, "low UMI"="grey"), guide="none")+
    scale_fill_manual("GABA", values=c("white","black"))+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  p = make_escheR(spe_sub) |> add_ground(var="precast_k9_1663_f", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="gad1", point_size = .5)
  p2 <- p+scale_color_manual("", values=c(cpList$earthy.pal1, "low UMI"="grey"), guide="none")+
    scale_fill_gradient("GAD1\nlog2\nCPM", low="white",high="black")+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  p = make_escheR(spe_sub) |> add_ground(var="precast_k9_1663_f", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="gad2", point_size = .5)
  p3 <- p+scale_color_manual("", values=c(cpList$earthy.pal1, "low UMI"="grey"), guide="none")+
    scale_fill_gradient("GAD2\nlog2\nCPM", low="white",high="black")+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  list(p1, p2, p3)
})
)
.nrow=6
.ncol=3
grobList = marrangeGrob(plist, layout_matrix=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T), top=NULL)
ggsave(file="plots/05_clustering/PRECAST/GABA-spot-plots_k9-1663.png", 
       grobList, height=12, width=9, bg="white")
cat("\nGABA spot plots saved to: plots/05_clustering/PRECAST/GABA-spot-plots_k9-1663.png\n")

#gaba violin plots
spe$gad1_more1 = spe$gad1>1
spe$gad2_more1 = spe$gad2>1

cdata = mutate(as.data.frame(colData(spe)), more1_gad1.gad2=paste(gad1_more1, gad2_more1),
               more1_gad1.gad2_label = paste0("GAD1>1 (", gad1_more1, ")\nGAD2>1 (", gad2_more1, ")"))

p1 <- ggplot(group_by(cdata, sample_id, precast_k9_1663_f, is_gaba, more1_gad1.gad2_label) %>% tally(), 
       aes(x=precast_k9_1663_f, y=n, color=is_gaba))+
  ggbeeswarm::geom_quasirandom(size=.5)+scale_color_manual(values=c("black","red3"), guide="none")+
  facet_wrap(vars(more1_gad1.gad2_label), scales="free_y", ncol=4)+
  labs(y="# spots per sample")+
  theme_bw()+theme(text=element_text(size=10), strip.text=element_text(size=12),
                   strip.background = element_blank(),
                   axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

### inhib subtype marker expression
spe_sub = spe[,spe$gad1>1 & spe$gad2>1]
cdata2 = as.data.frame(colData(spe_sub))
plot.genes = c("GAD1","GAD2","SLC6A1","SLC32A1",
               "LHX6","TAC1","PVALB","SST",
               "ADARB2","LAMP5","VIP","NPY")
for(i in plot.genes) {
  cdata2[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

cdata2 = group_by(cdata2, sample_id, precast_k9_1663_f, is_gaba) %>%
  summarise_at(all_of(plot.genes), mean) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="gene_name", values_to="avg_expr") %>%
  mutate(gene_name= factor(gene_name, levels=plot.genes))

p2 <- ggplot(cdata2, aes(x=precast_k9_1663_f, y=avg_expr, color=is_gaba))+
  ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c("black","red3"), guide="none")+
  facet_wrap(vars(gene_name), scales="free_y")+
  labs(title="GAD1>1 (TRUE) & GAD2>1 (TRUE)", y="per-sample avg. expr. (log2 CPM)")+
  theme_bw()+theme(text=element_text(size=10), strip.text=element_text(face="italic"),
                   strip.background=element_blank(),
                   axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank())

ggsave(file="plots/05_clustering/PRECAST/GABA-expr-violin_k9-1663.png", 
       grid.arrange(p1, p2, layout_matrix=rbind(c(1),c(2),c(2))),
       bg="white", width=12, height=12)
cat("\nGABA subtype expression plots saved to: plots/05_clustering/PRECAST/GABA-expr-violin_k9-1663.png\n")

#reduced dimension plots
### load in precast loadings
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1663", withDimnames=F) = reddim[spe$seurat_key,]

test_variables = c("sample_id","slide","sum_umi","seq","expr_chrM_ratio","sex","precast_k9_1663_f")

percVar.df = tibble::rownames_to_column(as.data.frame(getExplanatoryPCs(spe, dimred="PRECAST_1663", n_dimred = 15, 
		variables= test_variables)), 
	var="component") %>%
    tidyr::pivot_longer(all_of(test_variables), names_to="variable", values_to="percVariance")

percVar.df$component = as.numeric(sapply(strsplit(percVar.df$component, "_"), function(x) x[[2]]))

### spaghetti plot of % var per component
color.pal = RColorBrewer::brewer.pal(n=(length(test_variables)-1), "Set2")
p1 <- ggplot(mutate(percVar.df, variable=factor(variable, 
			levels=c("precast_k9_1663_f","sample_id","slide","seq","sex","sum_umi","expr_chrM_ratio"),
                        labels=c("cluster","sample","slide","seq","sex","sum UMI","mito. ratio"))), 
               aes(x=component, y=percVariance, color=variable))+
	geom_point()+geom_line(aes(group=variable))+
	ylim(0,100)+scale_x_continuous(breaks=c(1,5,10,15))+
	scale_color_manual(values=c("black",color.pal))+
	labs(title="PRECAST embeddings (k=9 n=1663): top variables explaining variance", 
		x="component", y="% of var explained", color="")+
	theme_minimal()+theme(text=element_text(size=10))

### silhouette plot
sil.results <- as.data.frame(approxSilhouette(reducedDim(spe, "PRECAST_1663"), clusters=spe$precast_k9_1663_f))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe$precast_k9_1663_f))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
	ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c(cpList$earthy.pal2, "low UMI"="grey"))+
	guides(color = guide_legend(override.aes = list(size = 1)))+
	labs(x="assigned PRECAST cluster", title="PRECAST embeddings (k=9 n=1663): silhouette")+
	theme_minimal()+theme(text=element_text(size=10))


### sample embeddings per cluster
med.df = as.data.frame(cbind.data.frame(colData(spe)[,c("sample_id","precast_k9_1663_f")], reducedDim(spe, "PRECAST_1663"))) %>%
	group_by(sample_id, precast_k9_1663_f) %>% summarise_all(median) %>% 
	tidyr::pivot_longer(all_of(paste("PRECAST",1:15, sep="_")), names_to="component", values_to="median_loading") %>%
	mutate(component=factor(component, levels=paste("PRECAST",1:15, sep="_")))

p3 <- ggplot(med.df, aes(x=precast_k9_1663_f, y=median_loading, color=precast_k9_1663_f))+
	ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c(cpList$earthy.pal2, "low UMI"="grey"))+
	facet_wrap(vars(component), ncol=4)+labs(y="per-sample median spot loading")+
	theme_minimal()+theme(text=element_text(size=8), panel.border=element_rect(fill=NA, color="black", linewidth=.3),
		axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1),
		legend.position="none", plot.margin=unit(c(4,8,4,4),"pt"))

ggsave(file="plots/05_clustering/PRECAST/reduced-dim-plots_k9-1663.png",
       grid.arrange(p1, p2, p3, layout_matrix=matrix(c(1,2,3,3))), 
	bg="white", width=8, height=12, units="in")
cat("\nReduced dim plots saved to: plots/05_clustering/PRECAST/reduced-dim-plots_k9-1663.png\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
