setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(edgeR)
	library(scuttle)
	library(scater)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")
low.res.pal = c("Astro"="#cfa45c","Micro.Vasc"="#911223",
                "Inhb"="#9377AC",
                "L2"="#5D9940","L3"="#5095CD",
                "L4"="#85A0A0",
                "L5"="#ddc94e","L6"="#E45C5F",
                "Oligo"="#D1C4B0")

load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")
sce_con = sce_pseudo
dim(sce_con)
sce_con$Age_death = as.numeric(as.character(sce_con$Age_death))

load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_bipolar_pseudo_indivID-low-res.Rdata")
sce_bd = sce_pseudo
dim(sce_bd)

sce_bd$seurat_low.res = factor(sce_bd$seurat_low.res, levels=c("Micro/Vasc", "Astro", "L2", "L3", "L4", "L5", "L6", "Oligo", "Inhb"),
	labels=c("Micro.Vasc", "Astro", "L2", "L3", "L4", "L5", "L6", "Oligo", "Inhb"))
sce_bd$Age_death = as.numeric(as.character(sce_bd$Age_death))


sce_pseudo = cbind(sce_con, sce_bd)


#filter by expression before recalculating norm counts
rowData(sce_con)$high_expr_group_individualID <- filterByExpr(sce_con, group = sce_con$individualID)
rowData(sce_con)$high_expr_group_cluster <- filterByExpr(sce_con, group = sce_con$seurat_low.res)

with(rowData(sce_con), table(high_expr_group_individualID, high_expr_group_cluster))

sce_keep.genes = rownames(sce_pseudo)[rowData(sce_pseudo)$high_expr_group_individualID & rowData(sce_pseudo)$high_expr_group_cluster]

# keep any genes in DE models
fdata = read.csv("processed-data/93_globus/MBv_n119_features.csv.gz", row.names=1)
de.genes = rownames(fdata)[fdata$DE_domain.sp|fdata$DE_domain.ct]

keep.genes = intersect(rownames(sce_pseudo), union(sce_keep.genes, de.genes))
length(keep.genes)


sce_pseudo = sce_pseudo[keep.genes,]
sce_pseudo <- addPerCellQC(sce_pseudo)

dxpal1 = cpList$dx.pal[c("NTC","BPD")]
names(dxpal1) = c("con","bd")
dxpal2 = cpList$dx.pal[c("NTC","BPD")]
names(dxpal2) = c("control","Bipolar Disorder")

cdata= as.data.frame(colData(sce_pseudo)) %>% 
        mutate(sex=factor(Biological_Sex, levels=c("female","male"), labels=c("F","M")),
                condition=factor(Disorder, levels=c("control","Bipolar Disorder"), labels=c("con","bd"))
        )
ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex, color=condition))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=dxpal1)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

sce_pseudo = sce_pseudo[,sce_pseudo$detected>5000]

tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=2)
min(x)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

m1 <- scran::modelGeneVar(sce_pseudo)
hvg.genes <- scran::getTopHVGs(m1, n=2000)

sce_pseudo <- runPCA(sce_pseudo, subset_row=hvg.genes, exprs_values="logcounts", name="PCA_HVG")
plotReducedDim(sce_pseudo, dimred="PCA_HVG", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)

geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
n1663.ids = rownames(sce_pseudo)[rowData(sce_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids) #1516
sce_pseudo <- runPCA(sce_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_SVG")
p1 <- plotReducedDim(sce_pseudo, dimred="PCA_SVG", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)
p2 <- plotReducedDim(sce_pseudo, dimred="PCA_SVG", ncomponents=2, colour_by = "Disorder", point_alpha=1)#+
        scale_color_manual("", values=dxpal2)

gridExtra::grid.arrange(p1, p2, ncol=2)
# issue is that the M/V pseudobulk for controls and BD are separated in PC space

# and it doesn't seem to reflect higher/lower vasc or microglia based on CLDN5 and CD74 expr
sce_pseudo$cldn5 = logcounts(sce_pseudo)[rowData(sce_pseudo)$gene_name=="CLDN5",]
sce_pseudo$cd74 = logcounts(sce_pseudo)[rowData(sce_pseudo)$gene_name=="CD74",]

plotReducedDim(sce_pseudo, dimred="PCA_SVG", ncomponents=2, colour_by = "cldn5", point_alpha=1)
