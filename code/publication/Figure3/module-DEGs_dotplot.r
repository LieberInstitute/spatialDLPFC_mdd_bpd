setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})
set.seed(123)


comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

cpList <- readRDS("plots/colorPalettes.rds")

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


refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")

col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)

# plot main figure modules PLP1
plot.genes = c("PLP1","MAG","ENPP2","TF","CLDN11","CNDP1",#"HSD11B1",
               "SGK1","SPP1",
               "LMNA","NEAT1","PGAM2","AQP1")
#               "WNK1","PTP4A2","ANP32B","MTURN")

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))
tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))


order1 = filter(refined.modules, TF=="PLP1", target %in% plot.genes[-1]) %>% 
  arrange(desc(importance)) %>% pull(target)
tmp$plot.genes = factor(tmp$gene_name, levels=c(rev(order1),"PLP1"))
tmp$gene_group = factor(tmp$gene_name, levels=plot.genes, 
                        labels=c(rep("A",6), rep("B",6)))
#, rep("B",4)))
#tmp$gene_group = factor(tmp$gene_group, levels=c("A","B","C"))




p1 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
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


tmp2 = filter(refined.modules, TF=="PLP1", target %in% plot.genes) %>%
  add_row(TF="PLP1", target="PLP1", importance=0) %>%
  mutate(y_lab=factor(target, levels=c(rev(order1),"PLP1")),
         gene_group = factor(target, levels=plot.genes, 
                             labels=c(rep("A",6), rep("B",6))))
#, rep("B",4))),
#         gene_group= factor(gene_group, levels=c("A","B","C")))


p1.2 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8), 
                        strip.background = element_blank(), strip.text = element_blank(),
                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))



gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% plot.genes]
names(gids) = rowData(sce_summ)[gids,"gene_name"]

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

tmp2 = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(gene_name= factor(gene_name, levels=c(rev(order1),"PLP1")),
        gene_group= factor(gene_name, levels=plot.genes, labels=c(rep("A",6),rep("B",6))))
#,rep("B",4))),
#	gene_group= factor(gene_group, levels=c("A","B","C")))

p1.1 = ggplot(tmp2, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        strip.background = element_blank(), strip.text = element_blank(),
                        legend.position="bottom", legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))

ggsave(file="plots/publication/Figure3/module-DEGs_dotplot_PLP1.pdf",
       marrangeGrob(grobs=list(p1, p1.1, p1.2), layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7), top=NULL),
       width=7, height=4)


# Vascular modules
plot.genes = c("IFITM3","C1R","SLCO4A1","TIMP1","OSMR","APOLD1","IL1R1",
#"CHI3L1","EDN1",
               "A2M","ABCG2","TNFSF10","RERGL","SLC38A5","RAMP2","MUSTN1")

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))
tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))


order2 = filter(refined.modules, TF=="IFITM3", target %in% plot.genes[2:7]) %>% 
  arrange(desc(importance)) %>% pull(target)
order2 = c("IFITM3", order2)

order3 = filter(refined.modules, TF=="A2M", target %in% plot.genes[9:14]) %>% 
  arrange(desc(importance)) %>% pull(target)
order3 = c("A2M", order3)

order4 = c(order2, order3)
tmp$plot.genes = factor(tmp$gene_name, levels=rev(order4))
tmp$gene_group = factor(tmp$gene_name, levels=plot.genes, 
                        labels=c(rep("A",7), rep("B",7)))


p2 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
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


tmp2 = bind_rows(filter(refined.modules, TF=="IFITM3", target %in% plot.genes) %>%
                   add_row(TF="IFITM3", target="IFITM3", importance=0),
                 filter(refined.modules, TF=="A2M", target %in% plot.genes) %>%
                   add_row(TF="A2M", target="A2M", importance=0)) %>%
  mutate(y_lab=factor(target, levels=rev(order4)),
         gene_group = factor(target, levels=plot.genes, 
                             labels=c(rep("A",7), rep("B",7))))


p2.2 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8), 
                        strip.background = element_blank(), strip.text = element_blank(),
                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))


gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% plot.genes]
names(gids) = rowData(sce_summ)[gids,"gene_name"]

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

tmp2 = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(gene_name= factor(gene_name, levels=rev(order4)),
        gene_group= factor(gene_name, levels=plot.genes, labels=c(rep("A",7),rep("B",7))))

p2.1 = ggplot(tmp2, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        strip.background = element_blank(), strip.text = element_blank(),
                        legend.position="bottom", legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))


ggsave(file="plots/publication/Figure3/module-DEGs_dotplot_IFITM3-A2M.pdf",
       arrangeGrob(grobs=list(p2, p2.1, p2.2), layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7), top=NULL),
       width=7, height=4)


# plot CD74
plot.genes = filter(refined.modules, TF=="CD74") %>% arrange(desc(importance)) %>% pull(target)
plot.genes = c("CD74", plot.genes)
# remove two genes that aren't LA or LR sig in seurat label results (since thats what is being plotted here)
plot.genes = setdiff(plot.genes, c("LYVE1","SLC2A5"))

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))

tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))

tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

p3 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_x_discrete(labels=c("L-A","M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2.5,2.5))+
  facet_grid(cols=vars(sex.group))+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  scale_size_continuous("adj. p",limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))

tmp2 = filter(refined.modules, TF=="CD74", target %in% plot.genes) %>%
  add_row(TF="CD74", target="CD74", importance=0) %>%
  mutate(y_lab=factor(target, levels=rev(plot.genes)))

p3.2 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8),
                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))


# snRNAseq ratio
gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% plot.genes]
names(gids) = rowData(sce_summ)[gids,"gene_name"]

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

tmp2 = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(gene_name= factor(gene_name, levels=levels(tmp$plot.genes)))

p3.1 = ggplot(tmp2, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        legend.position="bottom", legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))


ggsave(file="plots/publication/Figure3/module-DEGs_dotplot_CD74.pdf",
       arrangeGrob(grobs=list(p3, p3.1, p3.2), layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7)),
       width=7, height=6)

# inflamm module?
plot.genes = c("HSPA1A","HSPA1B","MAFF","CDKN1A","GADD45B",
	"MT1X","ANGPTL4","HILPDA","MT1A","ADM","CEBPD")
# GADD45B is in both modules and is pretty much equally as important is each module (so is ARID5A but I'm not highlighting that one for HSPA1A so nbd)
##  TF     target  importance regulation   rho     n
##  <chr>  <chr>        <dbl>      <int> <dbl> <int>
##1 MT1X   GADD45B       3.83          1 0.596     2
##2 HSPA1A GADD45B       3.74          1 0.561     2
##3 MT1X   ARID5A        1.27          1 0.643     2
##4 HSPA1A ARID5A        1.11          1 0.539     2
# deciding to only plot GADD45B with HSPA1A importance for simplicity
#plot.genes2 = filter(refined.modules, TF=="MT1X", target!="GADD45B") %>% arrange(desc(importance)) %>% pull(target)
#plot.genes = c(plot.genes, "MT1X", plot.genes2)

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5))) %>%
	mutate(logFC=round(logFC, 2))
tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))


tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
tmp$gene_group = factor(tmp$gene_name, levels=plot.genes,
                        labels=c(rep("A",5), rep("B",6)))


p4 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_x_discrete(labels=c("L-A","M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2.5,2.5))+
  facet_grid(cols=vars(sex.group), rows=vars(gene_group), scales="free_y", space="free_y")+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  scale_size_continuous("adj. p",limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8),
                   strip.background.y = element_blank(), strip.text.y = element_blank(),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))


tmp2 = bind_rows(filter(refined.modules, TF=="HSPA1A", target %in% plot.genes[1:5]) %>%
                   add_row(TF="HSPA1A", target="HSPA1A", importance=0),
                 filter(refined.modules, TF=="MT1X", target %in% plot.genes[6:11]) %>%
                   add_row(TF="MT1X", target="MT1X", importance=0)) %>%
  mutate(y_lab=factor(target, levels=rev(plot.genes)),
         gene_group = factor(target, levels=plot.genes,
                             labels=c(rep("A",5), rep("B",6))))


p4.2 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8),
                        strip.background = element_blank(), strip.text = element_blank(),
                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))

gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% plot.genes]
names(gids) = rowData(sce_summ)[gids,"gene_name"]

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

tmp2 = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(gene_name= factor(gene_name, levels=levels(tmp$plot.genes)),
	gene_group= factor(gene_name, levels=plot.genes, labels=c(rep("A",5),rep("B",6))))

p4.1 = ggplot(tmp2, aes(y=gene_name, x=mean.expr, fill=cellType))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+
  facet_grid(rows=vars(gene_group), scales="free_y", space="free_y")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
			strip.background = element_blank(), strip.text = element_blank(),
                        legend.position="bottom", legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))

ggsave(file="plots/publication/Figure3/module-DEGs_dotplot_HSPA1A-MT1X.pdf",
       marrangeGrob(grobs=list(p4, p4.1, p4.2), layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7)),
       width=7, height=4)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

