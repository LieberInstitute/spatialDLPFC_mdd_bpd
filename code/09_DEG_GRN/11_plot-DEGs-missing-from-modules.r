setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})

source("code/09_DEG_GRN/load_DEGs.r")
all.degs = unique(do.call(rbind, sigList)$gene_name)
cat("\nNumber of total DEGs:", length(all.degs), "\n")
#should be 503

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
missing.degs = setdiff(all.degs, union(refined.modules$TF, refined.modules$target))
cat("Number of DEGs missing from modules:", length(missing.degs), "\n")
#should be 63

# snRNAseq mean ratio bar plot
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")
gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% missing.degs]
names(gids) = rowData(sce_summ)[gids,"gene_name"]
cat("\nNumber of missing DEGs present in snRNAseq dataset:", length(gids), "\n")
#should be 62

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

bar.df = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), 
                           names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(cellType=factor(cellType, levels=c("Vasc","Astro","Oligo","Micro","InhN","ExcN")))

order1 = mutate(bar.df, is_nrn = cellType %in% c("InhN","ExcN")) %>% group_by(gene_name, is_nrn) %>%
  summarise(sum_mean= sum(mean.expr)) %>%
  tidyr::pivot_wider(names_from="is_nrn", values_from="sum_mean", names_prefix = "nrn_") %>%
  mutate(nrn.ratio= nrn_TRUE/nrn_FALSE) %>% arrange(nrn.ratio) %>%
  pull(gene_name)

# modify order by putting HB modules and ADAMTS1 module at end, followed by ones not present in snRNAseq dataset
#HB module
order2 = order1[c(grep("HB",order1), setdiff(1:length(order1), grep("HB",order1)))]
#ADAMTS1
adj = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")
c1 <- filter(adj, TF=="ADAMTS1", target %in% order1, importance>1)$target
order3 = c(c1, setdiff(order2, c1))
order3 = c(rev(c("TM4SF1","MGP","C11orf96","CRISPLD2","COL4A1","LMOD1","ADAMTS9")), setdiff(order2, c1))
#visual improvement on order
order4 = order3[c(1:12,14:15,16,13,17:53,55,56,58,57,60,59,54,61,62)]

#order3 = order2[c(c(4,10,13,15,16,17,29), setdiff(1:length(order2), c(4,10,13,15,16,17,29)))]
#order3[11:18] = c("ACY3","ITGAX","P2RY12","CHI3L2","GBP1","CD44","APLNR","SERTAD1")

order4 = c(setdiff(missing.degs, order4), order4)

bar.df$gene_name = factor(bar.df$gene_name, levels=order4)

col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)


p3.1 = ggplot(bar.df, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+scale_y_discrete(position="right")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.x=element_blank(), axis.text.y=element_text(size=6), axis.title.y=element_blank(),
    legend.position="bottom", legend.title=element_blank(),
    legend.key.size = unit(10,"pt"))

# plot dotplot
cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

## load in DE results and format
sigList <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  if(results_set=="smoothed-k9-1663") {
    comp_names= names(cpList$smoothed.bright)
    set_name = "sm"
  }
  if(results_set=="seurat-pc30") {
    comp_names= names(cpList$transfer.bright)[c(1:4,8,5:7)]
    set_name = "se"
  }
  
  la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
           adj.P.Val2=ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5),
           sex.group=factor(coef, levels=comparisons2),
           cluster="L-A", source=set_name) #%>%
    #filter(gene_id %in% la.degs$gene_id, adj.P.Val<.05)
  
  lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
           adj.P.Val2=ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5),
           cluster=factor(cluster, levels=comp_names), source=set_name,
           sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2)) #%>%
    #filter(gene_id %in% lr.degs$gene_id, adj.P.Val<.05)
  
  return(list("la"=lat, "lr"=lrt))
})

names(sigList) <- c("sm","se")
sigList = do.call(c, sigList)


sig.df = do.call(rbind, sigList) %>%
  mutate(cluster_source = paste(cluster, source))

all_clusters = c("L-A sm","L-A se","Micro.Vasc se","Astro se","L1 sm",
                 "L2 sm","L2.3 se","L3.4 sm","L4 se",
                 "Inhb se","L5 sm","L5 se","L6 sm","L6 se","WM sm","Oligo se")

sig.df = mutate(sig.df, cluster_source2 = factor(cluster_source, levels=all_clusters))

## filter to genes of interest
dot.df = filter(sig.df, gene_name %in% order4) %>% mutate(source=factor(source, levels=c("sm","se"))) 
dot.df$is_sig = dot.df$adj.P.Val2<.05
dot.df$gene_name = factor(dot.df$gene_name, levels=order4)

## format for overlap
dot.df$x_labels = factor(dot.df$cluster_source2, levels=all_clusters,
                      labels=c("L-A","L-A","Micro.Vasc","Astro","L1",
                               #"L2","L2.3","L3.4","L4",
                               "L2/3","L2/3","L3/4","L3/4",
                               "Inhb",
                               "L5","L5","L6","L6","WM/Oligo","WM/Oligo"))

dot.df$size2 = as.numeric(as.character(factor(paste(dot.df$source, dot.df$is_sig), levels=c("sm FALSE","se FALSE","sm TRUE","se TRUE"),
                   labels=c(1,1,4,3))))

p3 <- ggplot(dot.df, aes(x=x_labels, y=gene_name, fill=logFC, size=size2))+
  geom_count(aes(shape=source, color=is_sig))+
  scale_shape_manual(values=c("sm"=23, "se"=21))+
  #scale_size_manual("adj. p", values=c("FALSE"=1, "TRUE"=4))+ #for size=is_sig
  scale_size_identity("adj. p")+ #for size=size2
  scale_color_manual(values=c("FALSE"="grey", "TRUE"="black"))+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2.5,2.5))+
  facet_grid(cols=vars(sex.group), #rows=vars(gene_group), 
             scales="free_y", space="free_y")+
  scale_x_discrete(labels=c("L-A","M.V","Ast","L1","L2/3","L3/4","Inhb","L5","L6","WM/O"))+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8), axis.text.y=element_text(size=6),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))

ggsave(file="plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_degs-missing-from-modules.pdf",
	grid.arrange(p3, p3.1, layout_matrix=matrix(c(1,1,1,1,1,1,2), ncol=7)),
	height=8, width=8)
cat("\nMissing DEG plots saved to: plots/09_DEG_GRN/spe-n119_13162-no-lowUMI_degs-missing-from-modules.pdf\n")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
