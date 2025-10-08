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

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")
dim(spe)

cpList = readRDS("plots/colorPalettes.rds")
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)
stopifnot(identical(colnames(spe), rownames(cdata)))

spe$precast_k7_1626_f = factor(cdata$precast_k7_1626_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","missing"),
	labels=c("L1","L2","L3.4","L5","L6","WM","low UMI","missing"))
spe = spe[,spe$precast_k7_1626_f!="missing"]
dim(spe)

spe$precast_k7_1626_f = droplevels(spe$precast_k7_1626_f)
table(spe$precast_k7_1626_f, useNA="ifany")

#reduced dimension plots
### load in precast loadings
load("processed-data/05_clustering/PRECAST/srt_precast-conservative_k-7_n1626.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1626", withDimnames=F) = reddim[spe$seurat_key,]

test_variables = c("sample_id","slide","sum_umi","seq","expr_chrM_ratio","sex","precast_k7_1626_f")

percVar.df = tibble::rownames_to_column(as.data.frame(getExplanatoryPCs(spe, dimred="PRECAST_1626", n_dimred = 15, 
		variables= test_variables)), 
	var="component") %>%
    tidyr::pivot_longer(all_of(test_variables), names_to="variable", values_to="percVariance")

percVar.df$component = as.numeric(sapply(strsplit(percVar.df$component, "_"), function(x) x[[2]]))

### spaghetti plot of % var per component
color.pal = RColorBrewer::brewer.pal(n=(length(test_variables)-1), "Set2")
p1 <- ggplot(mutate(percVar.df, variable=factor(variable, 
			levels=c("precast_k7_1626_f","sample_id","slide","seq","sex","sum_umi","expr_chrM_ratio"),
                        labels=c("cluster","sample","slide","seq","sex","sum UMI","mito. ratio"))), 
               aes(x=component, y=percVariance, color=variable))+
	geom_point()+geom_line(aes(group=variable))+
	ylim(0,100)+scale_x_continuous(breaks=c(1,5,10,15))+
	scale_color_manual(values=c("black",color.pal))+
	labs(title="PRECAST embeddings (k=7 n=1626): top variables explaining variance", 
		x="component", y="% of var explained", color="")+
	theme_minimal()+theme(text=element_text(size=10))

### silhouette plot
sil.results <- as.data.frame(approxSilhouette(reducedDim(spe, "PRECAST_1626"), clusters=spe$precast_k7_1626_f))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(spe$precast_k7_1626_f))

p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
	ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c(cpList$smoothed.bright, "low UMI"="grey"))+
	guides(color = guide_legend(override.aes = list(size = 1)))+
	labs(x="assigned PRECAST cluster", title="PRECAST embeddings (k=7 n=1626): silhouette")+
	theme_minimal()+theme(text=element_text(size=10))


### sample embeddings per cluster
med.df = as.data.frame(cbind.data.frame(colData(spe)[,c("sample_id","precast_k7_1626_f")], reducedDim(spe, "PRECAST_1626"))) %>%
	group_by(sample_id, precast_k7_1626_f) %>% summarise_all(median) %>% 
	tidyr::pivot_longer(all_of(paste("PRECAST",1:15, sep="_")), names_to="component", values_to="median_loading") %>%
	mutate(component=factor(component, levels=paste("PRECAST",1:15, sep="_")))

p3 <- ggplot(med.df, aes(x=precast_k7_1626_f, y=median_loading, color=precast_k7_1626_f))+
	ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c(cpList$smoothed.bright, "low UMI"="grey"))+
	facet_wrap(vars(component), ncol=4)+labs(y="per-sample median spot loading")+
	theme_minimal()+theme(text=element_text(size=8), panel.border=element_rect(fill=NA, color="black", linewidth=.3),
		axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1),
		legend.position="none", plot.margin=unit(c(4,8,4,4),"pt"))

ggsave(file="plots/05_clustering/PRECAST/reduced-dim-plots_conservative_k7-1626.png",
       grid.arrange(p1, p2, p3, layout_matrix=matrix(c(1,2,3,3))), 
	bg="white", width=8, height=12, units="in")
cat("\nReduced dim plots saved to: plots/05_clustering/PRECAST/reduced-dim-plots_conservative_k7-1626.png\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
