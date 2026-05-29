setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(pheatmap)
	library(ggplot2)
	library(igraph)
	library(SpatialExperiment)
})

set.seed(123)

# load in GRN adjacency output for correlations
lg.mask = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")

# load in module sets
modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv")
modules2 = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

# read in full DEG corr mtx for all DEG module plotting
input.mtx = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv", row.names=1)

# pick out top predictor DEG network subset (top predictor network must have >6 components)
mod_subset = c("EEF1A1","EIF1","PLP1","CD74","COX4I1","GLUL","IFITM3","GRIN1",
	"SNHG14","CAMK2N1","GAD1","A2M","PRKAR1A","UQCRH","ADAMTS1","HSPA1A","MT1M","JUNB",
	"GFAP")

cat("\nNumber of INITIAL top predictor DEG network modules to highlight:", length(mod_subset), "\n")

each_module = unique(modules2$TF)
col_annot = data.frame("is_top"= as.character(each_module %in% mod_subset), row.names=each_module)
annot_colors= list("is_top"=c("FALSE"="white", "TRUE"="black"))
phm = pheatmap(input.mtx, clustering_method="ward.D2",
	breaks= seq(from = -1, to = 1, length.out = 101), 
	treeheight_col = 20, treeheight_row = 20, angle_col=90, fontsize=6,
	annotation_col=col_annot, annotation_row=col_annot, annotation_colors=annot_colors, 
	annotation_legend = FALSE, annotation_names_row = FALSE, annotation_names_col = FALSE,
	main="Pearson corr. of DEG module importance")


# now just the subsets I selected from top predictors
input.mtx = matrix(NA, nrow=length(mod_subset), ncol=length(mod_subset), dimnames = list(mod_subset, mod_subset))
for(i in 1:(nrow(input.mtx)-1)) {
  complete.corr = colnames(input.mtx)[(i+1):ncol(input.mtx)]
  i=rownames(input.mtx)[i]
  for(j in complete.corr) {
    tmp1 = filter(modules2, TF==i)
    tmp2 = filter(modules2, TF==j)
    geneList = setdiff(union(tmp1$target, tmp2$target),
                       c(i,j))
    #pull the raw adj (before any filtering) and format for correlation
    check = filter(lg.mask, TF %in% c(i,j), target %in% geneList) %>%
      mutate(modules_key=factor(TF, levels=c(i,j), labels=c("I","J"))) %>%
      select(modules_key, target, importance) %>% tidyr::pivot_wider(names_from="modules_key", values_from="importance", values_fill=0)
    corrij= cor(check$I, check$J)
    input.mtx[i,j] = corrij
    input.mtx[j,i] = corrij
  }
}

# plot top predictor modules heatmaps
phm1 = pheatmap(input.mtx, clustering_method="ward.D2",
                   breaks= seq(from = -1, to = 1, length.out = 101),
                   treeheight_col = 20, treeheight_row = 20, angle_col=90,
                main="Pearson corr. of module importance (DEG targets only)")

# identify highly correlated DEG modules and pick one
input.mtx[matrixStats::rowMaxs(input.mtx, na.rm=T)>.2,
          matrixStats::colMaxs(input.mtx, na.rm=T)>.2]
# keep UQCRH
cat("\nTwo top predictor networks represent highly correlated target gene modules...")
cat("\nRemove UQCRH and keep COX4I1...\n")
cat("\nRemove EIF1 and keep EE1F1A...\n")
exclude = c("UQCRH","EIF1")
mod_subset = setdiff(mod_subset, exclude)

cat("\nNumber of FINAL top predictor DEG network modules to highlight:", length(mod_subset), "\n")
cat("\nTop predictor DEG network modules to highlight:\n")
print(mod_subset)


# filter DEG modules to top predictor subset
tmp = filter(modules2, TF %in% mod_subset)

source("code/09_DEG_GRN/load_DEGs.r")
mbv.degs = unique(do.call(rbind, sigList)$gene_name)

cat("\nNumber of DEGs present in GRN output:", length(intersect(mbv.degs, lg.mask$target)),"\n")
cat("\nNumber of DEGs present in the", length(unique(modules$TF)), "interaction modules:", length(intersect(mbv.degs, modules$target)),"\n")
cat("\nNumber of DEGs present in the", length(unique(modules2$TF)), "DEG modules:", length(intersect(mbv.degs, modules2$target)),"\n")
cat("\nNumber of DEGs present in the", length(mod_subset), "top predictor DEG modules:", length(intersect(mbv.degs, tmp$target)),"\n\n")


# plot heatmap of expression for top predictor DEG modules
#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
spe_summ$seurat_label = factor(spe_summ$seurat_qual.genes_pc30.kweight50, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
                               labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))
seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")

sm.df = dotplotDF(spe_summ, mod_subset, swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(mod_subset)),
         clusters=factor(clusters, levels=levels(spe_summ$sample_id), labels=gsub(" ","\n", levels(spe_summ$sample_id))))
sm.df$condition = factor(substr(sm.df$clusters, start=0, stop=3), levels=c("NTC","MDD","BPD")) 
sm.df$sex = factor(substr(sm.df$clusters, start=5, stop=5), levels=c("F","M"))

#make labels
xlbs = levels(sm.df$clusters)
## keep only the odd numbered dx indicators
odd1 = seq(1, length(xlbs), by=2)
xlbs[-odd1] = substr(xlbs[-odd1], start=4, stop=11)

p1 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
  geom_count(aes(shape=sex, color=mean_expr_scaled, size=prop_spots))+
  scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=xlbs)+
  scale_size(range=c(1,5), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
         shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(logcounts,\nscaled)", size="Prop. of\nspots",
       y="top predictors", title="Top predictor DEG modules")+
  theme_minimal()+theme(axis.title.x=element_blank(), axis.text.x=element_text(size=7, hjust=0),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"))



pdf(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-subset_heatmap-dotplot.pdf", width=8, height=8)
plot(phm[[4]])
plot(phm1[[4]])
p1 #dotplot of expression
dev.off()
cat("\nSaved all plots to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_top-predictor-subset_heatmap-dotplot.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
