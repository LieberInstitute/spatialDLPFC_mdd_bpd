setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

cpList <- readRDS("plots/colorPalettes.rds")
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

# aucell violins ----
aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(aucell) = gsub("Regulon\\.for\\.","",colnames(aucell))

# now refined module subset (remove COX4I1 for redundancy and ADAMTS1 because too small)
mod_subset = c("A2M","IFITM3","CD74","HSPA1A","MT1X",
               "SNHG14","GLUL","CAMK2N1","GAD1","GRIN1","PRKAR1A",
               "UQCRH","APLP1","FTL","PLP1")

seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")
aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
                            labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$smoothed = factor(aucell$smoothed_k9_1663, levels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA"),
                         labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped"))

col.pal = cpList$transfer.bright[c(1:4,8,5:7)]
names(col.pal) = seurat_levels

aucell_long = tidyr::pivot_longer(aucell, all_of(mod_subset), names_to="module", values_to="AUCell") %>%
  mutate(module= factor(module, levels=mod_subset))


p1 <- ggplot(filter(aucell_long, module %in% c("PLP1","FTL","APLP1")),
             aes(x=seurat_label, y=AUCell, fill=seurat_label))+
  geom_violin(scale="width")+facet_grid(col=vars(module))+
  scale_fill_manual(values=col.pal, guide="none")+
  theme_minimal()+theme(axis.text.x=element_text(size=8))


p2 <- ggplot(filter(aucell_long, module %in% c("PLP1","FTL","APLP1"), smoothed!="dropped"),
             aes(x=smoothed, y=AUCell, fill=smoothed))+
  geom_violin(scale="width")+facet_grid(col=vars(module))+
  scale_fill_manual(values=cpList$smoothed.bright, guide="none")+
  theme_minimal()+theme(axis.text.x=element_text(size=8))

ggsave(file="plots/publication/Figure3/supp_Oligo-modules-AUCell_violin.pdf",
	grid.arrange(p1, p2, ncol=1),
	width=6, height=4)
# PLP1 module shown on volcanoes

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

results_set="seurat-pc30"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))


lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         sex.group=factor(coef, levels=comparisons2),
         cluster="L-A") 

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$transfer.bright)[c(1:4,8,5:7)]),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2))

lrt_oligo = filter(lrt, cluster %in% c("L5","L6","Oligo"))

df1 = bind_rows(mutate(lat, is_DEG= gene_name %in% la.degs$gene_name), 
                mutate(lrt_oligo, is_DEG= gene_name %in% lr.degs$gene_name)) %>% 
  mutate(cluster=factor(cluster, levels=c("L-A","L5","L6","Oligo"))) 

plp1.degs = c("PLP1", filter(refined.modules, TF=="PLP1")$target)
label_degs = c("PLP1","MAG","ENPP2","TF","CLDN11","CNDP1","SGK1","SPP1",
               "LMNA","NEAT1","PGAM2","WNK1","PTP4A2","ANP32B","MTURN")

plp1.volcano <- ggplot(df1, aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(color="grey", size=.5), dpi=300)+
  rasterize(geom_point(data=filter(df1, gene_name %in% plp1.degs, is_DEG==T), 
                       color="black", size=.5), dpi=300)+
  ggrepel::geom_text_repel(data=filter(df1, gene_name %in% label_degs, adj.P.Val<.05),
                           aes(label=gene_name), min.segment.length = 0, max.overlaps=Inf, size=2, color="purple")+
  xlim(-3.4,3.4)+scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  geom_hline(aes(yintercept=-log10(.05)), lty=2, color="red3")+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  theme_minimal()+theme(aspect.ratio=1, text=element_text(size=10))

ggsave(file="plots/publication/Figure3/supp_PLP1-module_volcano.pdf",
        plp1.volcano,
        width=6, height=6)




# FTL module DEG dotplot
plot.genes = c("RPL32","MAP4","RPS13","RPS27A",
               "TPT1","UBA52",
               "EIF1","DDIT4","ATF4","NFKBIA")

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))

tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))

order1 = filter(refined.modules, TF=="FTL", target %in% plot.genes) %>% 
  arrange(desc(importance)) %>% pull(target)
tmp$plot.genes = factor(tmp$gene_name, levels=rev(order1))
tmp$gene_group = factor(tmp$gene_name, levels=plot.genes, 
                        labels=c(rep("A",4), rep("B",2), rep("C",4)))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

p3 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_x_discrete(labels=c("L-A","M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2,2))+
  facet_grid(cols=vars(sex.group), rows=vars(gene_group), scales="free_y", space="free_y")+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  scale_size_continuous("adj. p",limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8),
                   strip.background.y = element_blank(), strip.text.y = element_blank(),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))


tmp2 = filter(refined.modules, TF=="FTL", target %in% plot.genes) %>%
  mutate(y_lab=factor(target, levels=rev(order1)),
         gene_group = factor(target, levels=plot.genes, 
                             labels=c(rep("A",4), rep("C",2), rep("B",4))))
p3.1 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8), 
                        strip.background = element_blank(), strip.text = element_blank(),
                        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))

ggsave(file="plots/publication/Figure3/supp_dotplot_DEG-GRNs_FTL.pdf",
       arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,2), ncol=7), top=NULL),
       width=6, height=4)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

