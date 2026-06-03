setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
  library(gridExtra)
})
set.seed(123)
setAutoBlockSize(1e9)

regulons <- read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons1 = filter(regulons, set_size>=10)$TF

#load aucell
aucell = read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)  

# plot avg aucell for each regulon for each sample
df1 = group_by(aucell, condition, sex, sample_id) %>%
  summarise_at(regulons1, mean) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")),
         sex=factor(sex, levels=c("F","M")))
  
df1 = tidyr::pivot_longer(df1, all_of(regulons1), names_to="regulon", values_to="avg_aucell")

reg_order = group_by(df1, regulon) %>% summarise(med_aucell=median(avg_aucell)) %>% 
  arrange(med_aucell) %>% pull(regulon)

p1 <- ggplot(mutate(df1, regulon = factor(regulon, levels=reg_order)), 
             aes(x=regulon, y=avg_aucell))+
  ggbeeswarm::geom_quasirandom(size=.5)+
  scale_y_continuous(labels=function(x) sprintf("%.2f", x))+
  labs(x="", y="avg. AUCell (per donor)", title="Initial regulons (n=119)")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


p2 <- ggplot(mutate(df1, regulon = factor(regulon, levels=reg_order)), 
       aes(x=regulon, y=avg_aucell, group=sample_id))+
  geom_line()+scale_y_continuous(labels=function(x) sprintf("%.2f", x))+
  labs(x="", y="avg. AUCell (per donor)", title="Initial regulons (n=119)")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))



# look at 3 AP1 regulons
ap1.df = filter(ungroup(df1), regulon %in% c("FOS","JUNB","FOSB")) %>% 
  group_by(regulon) %>% mutate(rank_reg=rank(-avg_aucell)) %>%
  group_by(sample_id) %>% mutate(med_rank=median(rank_reg))

ap1.df$y_labs = factor(ap1.df$sample_id, levels=arrange(ap1.df, med_rank) %>% pull(sample_id) %>% unique())

p3 <- ggplot(ap1.df, aes(x=y_labs, y=avg_aucell, color=regulon))+
  geom_point()+geom_vline(aes(xintercept=10.5))+
  labs(x="sample_id (ordered by AP1 regulons)", y="avg. AUCell", color="Initial\nregulon",
       title="AP1 donor-enriched regulons")+
  theme_bw()+theme(#axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
#                     legend.position="inside", legend.position.inside = c(.8,.7),
		   legend.position="bottom", axis.text.x=element_blank(),
                   panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# use top3 targets to identify npas4 (cytoscape for plot)
# plot npas4 prop detected
load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot-no-lowUMI_sample-id.Rdata")
cpList <- readRDS("plots/colorPalettes.rds")

npas4.df = data.frame(sample_id=spe_summ$sample_id, 
                  npas4_prop.detected=assay(spe_summ, "logcounts.prop.detected")[rowData(spe_summ)$gene_name=="NPAS4",]) %>%
  mutate(y_labs=factor(sample_id, levels=arrange(ap1.df, med_rank) %>% pull(sample_id) %>% unique())) %>%
  left_join(distinct(aucell[,c("sample_id","condition","sex")]))

p4 <- ggplot(npas4.df, aes(x=y_labs, y=npas4_prop.detected, color=condition, shape=sex))+
  geom_point()+scale_color_manual(values=cpList$dx.pal)+
  scale_shape_manual(values=c("F"=19, "M"=17))+
  ylim(0,1)+geom_vline(aes(xintercept=10.5))+
  labs(x="sample_id (ordered by AP1 regulons)", y="prop. spots with NPAS4", 
       title="NPAS4 outliers in donor-enriched AP1 regulons")+
  theme_bw()+theme(#axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
#        legend.position="inside", legend.position.inside = c(.8,.6), 
        legend.position="bottom", axis.text.x=element_blank(),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())



# pull top 16 to plot 10 outliers (removed) and give context for the rest
plot.npas4.outliers = levels(npas4.df$y_labs)[1:16]


#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots
cdata2 = cdata[cdata$smoothed_k9_1663_f!="low UMI",]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","drop","drop"))

spe[["NPAS4"]] = logcounts(spe)[rowData(spe)$gene_name=="NPAS4",]


max.val = max(colData(spe)[["NPAS4"]])
max.valr = round(max.val,1)
if(max.valr<max.val) max.valr=max.valr+.1

plist1 <- lapply(plot.npas4.outliers, function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="NPAS4", point_size = .5)
  p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"), guide="none")+
    scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
    labs(subtitle=unique(spe_sub$sample_id))+
    theme(plot.subtitle=element_text(size=9))
  return(rasterize(p2, dpi=200))
})

plist2 <- lapply(plot.npas4.outliers[1:4], function(x) {
  spe_sub = spe[,spe$sample_id==x]
  p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                         stroke=.3, point_size = .5) |> 
    add_fill(var="NPAS4", point_size = .5)
  p2 <- p+scale_color_manual(values=c(cpList$smoothed.light, "drop"="grey"))+
    scale_fill_gradient(limits=c(0,max.valr), low="white",high="black")+
    labs(subtitle=unique(spe_sub$sample_id))+
    theme(plot.subtitle=element_text(size=9))
  return(rasterize(p2, dpi=200))
})

.nrow=4
.ncol=4
lm=matrix(seq_len(.nrow*.ncol), nrow = .nrow, ncol = .ncol, byrow = T)

pdf(file="plots/10_SCENIC/spe-n119_top-16_NPAS4-expr_spot-plot.pdf", height=8, width=9)
plot(arrangeGrob(grobs=plist1, layout_matrix=lm, top="NPAS4 expression (fill scale fixed, can compare across sections)"))
plot(arrangeGrob(grobs=plist2, ncol=2, top="NPAS4 expression with example legends"))
dev.off()
cat("\nNPAS4 expression spot plots saved to: plots/10_SCENIC/spe-n119_top-16_NPAS4-expr_spot-plot.pdf\n")


# now plot true regulon results (after removing 10 npas4 outliers)
aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)

regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons2 = filter(regulons, set_size>=10)$TF


df2 = group_by(aucell, condition, sex, sample_id) %>%
  summarise_at(regulons2, mean) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")),
         sex=factor(sex, levels=c("F","M")))

df2 = tidyr::pivot_longer(df2, all_of(regulons2), names_to="regulon", values_to="avg_aucell")

reg_order = group_by(df2, regulon) %>% summarise(med_aucell=median(avg_aucell)) %>% 
  arrange(med_aucell) %>% pull(regulon)

p5 <- ggplot(mutate(df2, regulon = factor(regulon, levels=reg_order)), 
             aes(x=regulon, y=avg_aucell))+
  ggbeeswarm::geom_quasirandom(size=.5)+
  labs(x="", y="avg. AUCell (per donor)", title="Final regulons (n=109)")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


p6 <- ggplot(mutate(df2, regulon = factor(regulon, levels=reg_order)), 
             aes(x=regulon, y=avg_aucell, group=sample_id))+
  geom_line()+
  labs(x="", y="avg. AUCell (per donor)", title="Final regulons (n=109)")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


ggsave(file="plots/10_SCENIC/spe-n119_identifying-NPAS4-outlier-donors.pdf",
       marrangeGrob(grobs=list(p1, p2, p3, p4, p5, p6), ncol=1, nrow=2, top=NULL),
       height=6, width=6)
cat("\nProcess plots saved to: plots/10_SCENIC/spe-n119_identifying-NPAS4-outlier-donors.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
