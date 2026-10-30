setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(edgeR)
	library(scuttle)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(bluster)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
names(cpList$dx.pal)[3] = "BD"
low.res.pal = c("Ast"="#cfa45c","M.V"="#911223",
                "Inb"="#9377AC",
                "L2"="#5D9940","L3"="#5095CD",
                "L4"="#85A0A0",
                "L5"="#ddc94e","L6"="#E45C5F",
                "Olg"="#D1C4B0")

load("processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res.Rdata")
dim(sce_pseudo)
rowData(sce_pseudo)$gene_name = rowData(sce_pseudo)$featurekey

# remove NA columns after pseudobulking
colData(sce_pseudo)$Channel <- NULL
colData(sce_pseudo)$key <- NULL

# remove the many donor-level QC metadata (for simplicity)
colData(sce_pseudo) <- colData(sce_pseudo)[,c(1:2,7,12:16)]
#colData(sce_pseudo) <- colData(sce_pseudo)[,c(1:3,8,13:17)] # for when/if including channel
colData(sce_pseudo)$individualID = droplevels(sce_pseudo$individualID)

# filter genes
rowData(sce_pseudo)$high_expr_group_indivID <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$seurat_low.res)

with(rowData(sce_pseudo), table(high_expr_group_indivID, high_expr_group_cluster))
sce_keep.genes = rownames(sce_pseudo)[rowData(sce_pseudo)$high_expr_group_indivID & rowData(sce_pseudo)$high_expr_group_cluster]
cat("\nGenes to keep from edgeR filtered gene set:", length(sce_keep.genes), "\n")

#sce_pseudo = sce_pseudo[sce_keep.genes,]

sce_save = sce_pseudo
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")
rowData(sce_pseudo)$high_expr_group_individualID <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$seurat_low.res)
sce_keep.genes2 = rownames(sce_pseudo)[rowData(sce_pseudo)$high_expr_group_individualID & rowData(sce_pseudo)$high_expr_group_cluster]

keep.genes = union(sce_keep.genes, sce_keep.genes2)
cat("\nNumber of genes to keep (union of SZBD control edgeR filters and Maitra edgeR filters):", length(keep.genes), "\n")
sce_pseudo = sce_save[keep.genes,]

#fdata = read.csv("processed-data/93_globus/MBv_n119_features.csv.gz", row.names=1)
#de.genes = rownames(fdata)[fdata$DE_domain.sp | fdata$DE_domain.ct]
#cat("Genes to keep from MBv SRT (used as DE input in either annot):", length(de.genes), "\n")

#keep.genes = union(sce_keep.genes, de.genes)
#keep.genes = intersect(rownames(sce_pseudo), keep.genes)
#cat("\nNumber of genes to keep (union of edgeR filters and MBv DE input):", length(keep.genes), "\n")

## this is too many genes to keep so instead I will just do the same genes from Maitra
#cat("\n!!!! Don't do that, instead just keep the same genes that were kept in the Maitra dataset !!!\n")
#sce_save = sce_pseudo
#load("processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res_norm-filt.Rdata")
#maitra_keep.genes = rownames(sce_pseudo)
#keep.genes = intersect(maitra_keep.genes, rownames(sce_save))
#length(keep.genes)
#sce_pseudo <- sce_save[keep.genes,]


ribo.genes = rownames(sce_pseudo)[grep("RPS|RPL",rowData(sce_pseudo)$gene_name)]
cat("\nNumber of ribo genes:\n")
length(ribo.genes)


sce_pseudo <- addPerCellQC(sce_pseudo, subsets = list(ribo=ribo.genes))


# plot QC
cdata= as.data.frame(colData(sce_pseudo))

p1 <- ggplot(cdata, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_wrap(vars(seurat_low.res), ncol=9, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of nuclei per pseudobulked sample", y="nspots")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())


p2 <- ggplot(cdata, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata, aes(x=condition, y=subsets_ribo_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+ylim(0,20)+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of RPS|RPL reads per pseudobulked sample", y="subsets_ribo_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_unfiltered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_unfiltered_QC-metrics.png\n")

# normalize
tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=10)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

# heatmap of layer markers
source("code/06_pseudobulk/dlpfc_genes.r")
hm1 <- plotHeatmap(sce_pseudo, features=unlist(dlpfc.genes), swap_rownames="gene_name",
	cluster_rows=F, colour_columns_by=c("seurat_low.res","detected"), column_annotation_colors=list("seurat_low.res"=low.res.pal),
	scale=T, center=T, annotation_legend=F, legend=F, fontsize=6, treeheight_col=15)

#PCA before filtering reveals one component dominated by low detected genes samples
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
n1663.ids = rownames(sce_pseudo)[rowData(sce_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids)

sce_pseudo <- runPCA(sce_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1663")

# PCA with HVGs looks better and still has batch/sample element
m1 <- scran::modelGeneVar(sce_pseudo)
hvg.genes <- scran::getTopHVGs(m1, n=2000)
rowData(sce_pseudo)$HVGs = rownames(sce_pseudo) %in% hvg.genes

sce_pseudo <- runPCA(sce_pseudo, subset_row=hvg.genes, exprs_values="logcounts", name="PCA_HVG")

exp.vars = c("seurat_low.res","individualID",
             "batch",
             "condition","sex",
             "Age_death","PMI",#"RIN",
	"subsets_ribo_percent",
             "sum","detected","ncells")

exp.vars.colors = c("#FB8072", "#80B1D3",
                    "#BC80BD", #"#B3DE69","#FFFFB3",
                    "#8DD3C7", "#FDB462",
                    "#FCCDE5", "#BEBADA", "#FFFFB3",
                    "#CCEBC5", "#D9D9D9","black")
names(exp.vars.colors) = exp.vars

p2 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_HVG", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(sce_pseudo, dimred="PCA_HVG", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)
p4 <- plotPCA(sce_pseudo, dimred="PCA_HVG", ncomponents=4, colour_by = "seurat_low.res")+
        scale_color_manual("", values=low.res.pal)

#retroactively modify linewidth and point size
change_linewidth = function(.plot) {
        m = ggplot_build(.plot)
        m$data[[1]]$size = 2
        m$data[[2]]$linewidth = 1
        m1 = ggplot_gtable(m)
        return(m1)
}
#q1 <- change_linewidth(p1)
q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

#p5 <- plotPCA(sce_pseudo, dimred="PCA_HVG", ncomponents=4, colour_by = "detected")+
#        scale_color_viridis_c("detected\ngenes")
#q5 = ggplot_build(p5)
#q5$data[[2]]$size = 1
#q5 = ggplot_gtable(q5)

bp1 <- ggplot(group_by(cdata, condition, sex, seurat_low.res) %>% summarise(n_total=sum(ncells)) %>% mutate(cond_sex=factor(paste(condition, sex), 
	levels=c("NTC F","BD F","NTC M","BD M"))), 
       aes(x=cond_sex, y=n_total, fill=seurat_low.res))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("# nuclei", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("", values=low.res.pal)+
  labs(title="Cell type abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(),
	legend.key.size=unit(8,"pt"))


# silhouette
sil.results <- as.data.frame(approxSilhouette(reducedDim(sce_pseudo, "PCA_HVG"), clusters=sce_pseudo$seurat_low.res))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(sce_pseudo$seurat_low.res))
p5 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
	ggbeeswarm::geom_quasirandom()+
	scale_color_manual(values=low.res.pal)+
	theme_minimal()

ggsave("plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_unfiltered_PCA-HVG-eval.png",
	gridExtra::grid.arrange(hm1[[4]], q2, p3, q4, bp1, p5, layout_matrix=cbind(c(1,3,5),c(2,4,6))),
        bg="white", height=10, width=12, units="in")
cat("\nPCA eval plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_unfiltered_PCA-HVG-eval.png\n")


save(sce_pseudo, file="processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm.Rdata")
cat("\nNormalized pseudobulk sce saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm.Rdata\n")

##### stop early to evaluate filters then run again
#stop("Evaluate plots then continue")
#####

#apply sample filter
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), levels=c("NTC F","BD F","NTC M","BD M"))
cat("\nSample number per group before QC filters:\n")
table(cdata[,c("seurat_low.res","cond_sex")])

cdata2 = cdata[cdata$ncells>=5,]
cdata2 = cdata2[cdata2$detected>5000,]

### extra filters
#rownames(cdata2) = paste(cdata2$individualID, cdata2$Channel, cdata2$seurat_low.res, sep="_")
#remove.nrn = rownames(cdata2)[cdata2$seurat_low.res %in% c("L2","L3","L4","L5","L6","Inb") &
#	cdata2$sum<1e5]
#cdata2 = cdata2[!rownames(cdata2) %in% remove.nrn,]

#cdata2 = cdata2[cdata2$individualID!="BD1",]
### based on silhouette score after the above first round of filters, there are several M/V samples that don't fit well
### these happen to be all MDD M M.V samples, which when removed equalize the number of case/control M samples
#rownames(cdata2) = paste(cdata2$Sample, cdata2$seurat_low.res, sep="_")
#remove.mv = c("M28_M.V","M34_M.V","M26_M.V","M14_M.V")
#cdata2 = cdata2[!rownames(cdata2) %in% remove.mv,]
### additional detected filter for astro and oligo
#remove.glia = rownames(cdata2)[cdata2$seurat_low.res %in% c("Ast","Olg") & cdata2$detected<8000]
#cdata2 = cdata2[!rownames(cdata2) %in% remove.glia,]
## last additional detected filter for L6
#remove.l6 = rownames(cdata2)[cdata2$seurat_low.res=="L6" & cdata2$detected<10000]
#cdata2 = cdata2[!rownames(cdata2) %in% remove.l6,]

cat("\nSample number per group after QC filters:\n")
table(cdata2[,c("seurat_low.res","cond_sex")])

# now filter sce object
sce_pseudo = sce_pseudo[,sce_pseudo$ncells>=5]
sce_pseudo = sce_pseudo[,sce_pseudo$detected>5000]
#sce_pseudo = sce_pseudo[,sce_pseudo$individualID!="BD1"]
#colnames(sce_pseudo) = paste(sce_pseudo$Sample, sce_pseudo$seurat_low.res, sep="_")
#colnames(sce_pseudo) = paste(sce_pseudo$individualID, sce_pseudo$Channel, sce_pseudo$seurat_low.res, sep="_")
#sce_pseudo = sce_pseudo[,rownames(cdata2)]

dim(sce_pseudo)
#sce_pseudo$individualID = droplevels(sce_pseudo$individualID)

#renorm after removing samples
tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=10)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

#boxplots of QC metrics by condition and sex
p1 <- ggplot(cdata2, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+                     
  facet_wrap(vars(seurat_low.res), ncol=9, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of nuclei per pseudobulked sample", y="ncells")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank(),
	axis.text.y=element_text(size=8))

p2 <- ggplot(cdata2, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata2, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata2, aes(x=condition, y=subsets_ribo_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+ylim(0,20)+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of ribo reads per pseudobulked sample", y="subsets_ribo_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_QC-metrics.png\n")

#heatmap of layer markers
hm2 <- plotHeatmap(sce_pseudo, features=unlist(dlpfc.genes), swap_rownames="gene_name",
        cluster_rows=F, colour_columns_by=c("seurat_low.res","detected"), column_annotation_colors=list("seurat_low.res"=low.res.pal),
        scale=T, center=T, annotation_legend=F, legend=F, fontsize=6, treeheight_col=15)

# PCA with HVGs looks better and still has batch/sample element
m1 <- scran::modelGeneVar(sce_pseudo)
hvg.genes <- scran::getTopHVGs(m1, n=2000)

sce_pseudo <- runPCA(sce_pseudo, subset_row=hvg.genes, exprs_values="logcounts", name="PCA_HVG")

p2 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_HVG", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(sce_pseudo, dimred="PCA_HVG", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)
p4 <- plotPCA(sce_pseudo, dimred="PCA_HVG", ncomponents=4, colour_by = "seurat_low.res")+
        scale_color_manual("", values=low.res.pal)

q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

#p5 <- plotPCA(sce_pseudo, dimred="PCA_HVG", ncomponents=4, colour_by = "detected")+
#        scale_color_viridis_c("detected\ngenes")
#q5 = ggplot_build(p5)
#q5$data[[2]]$size = 1
#q5 = ggplot_gtable(q5)

bp1 <- ggplot(group_by(cdata2, cond_sex, seurat_low.res) %>% summarise(n_total=sum(ncells)),
       aes(x=cond_sex, y=n_total, fill=seurat_low.res))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("# nuclei", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("", values=low.res.pal)+
  labs(title="Cell type abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(),
        legend.key.size=unit(8,"pt"))

# silhouette
sil.results <- as.data.frame(approxSilhouette(reducedDim(sce_pseudo, "PCA_HVG"), clusters=sce_pseudo$seurat_low.res))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(sce_pseudo$seurat_low.res))
p5 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
        ggbeeswarm::geom_quasirandom()+
        scale_color_manual(values=low.res.pal)+
        theme_minimal()


ggsave("plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_PCA-HVG-eval.png",
        gridExtra::grid.arrange(hm2[[4]], q2, p3, q4, bp1, p5, layout_matrix=cbind(c(1,3,5),c(2,4,6))),
        bg="white", height=10, width=12, units="in")
cat("\nPCA eval plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_PCA-HVG-eval.png\n")


# save with SVG PCA too
set.seed(123)
sce_pseudo <- runPCA(sce_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1663")

p2 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_1663", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(sce_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)
p4 <- plotPCA(sce_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "seurat_low.res")+
        scale_color_manual("", values=low.res.pal)

q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

sil.results <- as.data.frame(approxSilhouette(reducedDim(sce_pseudo, "PCA_1663"), clusters=sce_pseudo$seurat_low.res))
sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
sil.results$closest <- factor(sil.results$closest, levels=levels(sce_pseudo$seurat_low.res))
p5 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
        ggbeeswarm::geom_quasirandom()+
        scale_color_manual(values=low.res.pal)+
        theme_minimal()

ggsave("plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_PCA-SVG-eval.png",
        gridExtra::grid.arrange(hm2[[4]], q2, p3, q4, bp1, p5, layout_matrix=cbind(c(1,3,5),c(2,4,6))),
        bg="white", height=10, width=12, units="in")
cat("\nPCA eval plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/SZBD-low-res_filtered_PCA-SVG-eval.png\n")

save(sce_pseudo, file="processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm-filt.Rdata")
cat("\nFiltered, normalized pseudobulk sce saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm-filt.Rdata\n")


cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

