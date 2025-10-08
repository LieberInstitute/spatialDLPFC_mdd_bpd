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

seurat_pc = "pc20"
#low.res.pal = c("Micro/Vasc"="#911223","Astro"="#cfa45c",
#        "L2"="#5D9940", "L3"="#5095CD",
#        "L4"="#85A0A0",
#        "L5"="#ddc94e","L6"="#E45C5F",
#        "Oligo"="#D1C4B0",
#        "Inhb"="#9377AC")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")

res = read.csv(paste0("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-",
        seurat_pc, "_red-precast-kweight-50-low-res.csv"), row.names=1)
#stopifnot(identical(colnames(spe), rownames(res)))

# modify label for low UMI GABA and L4
# remove low UMI cluster for L4 and Inhb
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)
cdata = cdata[rownames(res),]
cdata$predicted.id = res$predicted.id
remove.spots = cdata$smoothed_k7_1626=="low UMI" & cdata$predicted.id %in% c("L4","Inhb")

cdata[remove.spots & cdata$predicted.id=="L4","predicted.id"] = "L4 (low UMI)"
cdata[remove.spots & cdata$predicted.id=="Inhb","predicted.id"] = "Inhb (low UMI)" 

spe = spe[,rownames(cdata)]
#res = res[colnames(spe),]

spe$seurat_label = factor(cdata$predicted.id, levels=c("Micro.Vasc","Astro","L2","L3","L4","L4 (low UMI)","L5","L6","Oligo","Inhb (low UMI)","Inhb"),
        ########## merging pc30 L2 and L3
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L4 (low UMI)","L5","L6","Oligo","Inhb (low UMI)","Inhb"))
        #labels=c("Micro/Vasc","Astro","L2/3","L2/3","L4","L5","L6","Oligo","Inhb"))
        ##########
table(spe$seurat_label)

######### merging pc30 L2 and L3
#seurat_pc = "pc30-L2-L3-merge"


cpList = readRDS("plots/colorPalettes.rds")
#cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(colnames(spe), rownames(cdata)))
spe$smoothed_k7_1626 = factor(cdata$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
	labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))
#table(spe$precast_k9_1663_f)

#check GAD1 and GAD2
spe$is_gaba = factor(ifelse(spe$seurat_label %in% c("Inhb","Inhb (low UMI)"), as.character(spe$seurat_label), "other"),
	levels=c("other","Inhb (low UMI)","Inhb"))
spe$gad1 = logcounts(spe)[rowData(spe)$gene_name=="GAD1",]
spe$gad2 = logcounts(spe)[rowData(spe)$gene_name=="GAD2",]

#gaba spot plots
sub_samples = paste0("V13B23-", c("329_A1","332_A1","309_D1","352_B1","352_D1","308_A1"))

plist <- do.call(c, lapply(sub_samples, function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k7_1626", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="is_gaba", point_size = .5)
  p1 <- p+scale_color_manual("", values=c(cpList$smoothed.light, "low UMI"="grey"), guide="none")+
    scale_fill_manual("GABA", values=c("other"="white","Inhb (low UMI)"="red","Inhb"="black"))+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k7_1626", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="gad1", point_size = .5)
  p2 <- p+scale_color_manual("", values=c(cpList$smoothed.light, "low UMI"="grey"), guide="none")+
    scale_fill_gradient("GAD1\nlog2\nCPM", low="white",high="black")+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k7_1626", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="gad2", point_size = .5)
  p3 <- p+scale_color_manual("", values=c(cpList$smoothed.light, "low UMI"="grey"), guide="none")+
    scale_fill_gradient("GAD2\nlog2\nCPM", low="white",high="black")+
    theme(text=element_text(size=10), legend.key.size = unit(8,"pt"))
  list(p1, p2, p3)
})
)
.nrow=6
.ncol=3
grobList = marrangeGrob(plist, layout_matrix=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T), top=NULL)
ggsave(file=paste0("plots/05_clustering/Seurat/GABA-spot-plots_MBv-conservative_transfer-", seurat_pc, ".png"), 
       grobList, height=12, width=9, bg="white")
cat("\nGABA spot plots saved to:",paste0("plots/05_clustering/Seurat/GABA-spot-plots_MBv-conservative_transfer-", seurat_pc, ".png"),"\n")

#gaba violin plots
spe$gad1_more1 = spe$gad1>1
spe$gad2_more1 = spe$gad2>1

cdata = mutate(as.data.frame(colData(spe)), more1_gad1.gad2=paste(gad1_more1, gad2_more1),
               more1_gad1.gad2_label = paste0("GAD1>1 (", gad1_more1, ")\nGAD2>1 (", gad2_more1, ")"))

p1 <- ggplot(group_by(cdata, sample_id, seurat_label, is_gaba, more1_gad1.gad2_label) %>% tally(), 
       aes(x=seurat_label, y=n, color=is_gaba))+
  ggbeeswarm::geom_quasirandom(size=.5)+scale_color_manual(values=c("grey","red3","black"), guide="none")+
  facet_wrap(vars(more1_gad1.gad2_label), scales="free_y", ncol=4)+
  labs(y="# spots per sample")+
  theme_bw()+theme(text=element_text(size=10), strip.text=element_text(size=12),
                   strip.background = element_blank(),
                   axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

### inhib subtype marker expression
#both.false = !spe$gad1>1 & !spe$gad2>1
#spe_sub = spe[,!both.false]
spe_sub = spe[,spe$gad1_more1 & spe$gad2_more1]
cdata2 = as.data.frame(colData(spe_sub))
plot.genes = c("GAD1","GAD2","SLC6A1","SLC32A1",
               "LHX6","TAC1","PVALB","SST",
               "ADARB2","LAMP5","VIP","NPY")
for(i in plot.genes) {
  cdata2[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

cdata2 = group_by(cdata2, sample_id, seurat_label, is_gaba) %>%
  summarise_at(all_of(plot.genes), mean) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="gene_name", values_to="avg_expr") %>%
  mutate(gene_name= factor(gene_name, levels=plot.genes))

p2 <- ggplot(cdata2, aes(x=seurat_label, y=avg_expr, color=is_gaba))+
  ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=c("grey","red3","black"), guide="none")+
  facet_wrap(vars(gene_name), scales="free_y")+
  labs(title="GAD1>1 (TRUE) & GAD2>1 (TRUE)", y="per-sample avg. expr. (log2 CPM)")+
  theme_bw()+theme(text=element_text(size=10), strip.text=element_text(face="italic"),
                   strip.background=element_blank(),
                   axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank())

ggsave(file=paste0("plots/05_clustering/Seurat/GABA-expr-violin_MBv-conservative_transfer-", seurat_pc, ".png"), 
       grid.arrange(p1, p2, layout_matrix=rbind(c(1),c(2),c(2))),
       bg="white", width=12, height=12)
cat("\nGABA subtype expression plots saved to:", paste0("plots/05_clustering/Seurat/GABA-expr-violin_MBv-conservative_transfer-", seurat_pc,".png"),"\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
