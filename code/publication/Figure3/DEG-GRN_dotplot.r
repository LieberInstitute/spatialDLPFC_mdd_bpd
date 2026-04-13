setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
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

# plot main figure modules PLP1, IFITM3, A2M
plot.genes = c("PLP1","MAG","ENPP2","TF","CLDN11","CNDP1","SGK1","SPP1",
		"LMNA","HSD11B1","NEAT1","PGAM2","AQP1",
	"IFITM3","C1R","SLCO4A1","TIMP1","OSMR","APOLD1","CHI3L1","EDN1","IL1R1",
	"A2M","ABCG2","TNFSF10","RERGL","SLC38A5","RAMP2","MUSTN1")

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))

tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))

tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

p1 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_x_discrete(labels=c("L-A","M.V","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2,2))+
  facet_grid(cols=vars(sex.group))+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  scale_size_continuous("adj. p",limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))

tmp2 = bind_rows(filter(refined.modules, TF=="PLP1", target %in% plot.genes) %>%
  add_row(TF="PLP1", target="PLP1", importance=0),
	filter(refined.modules, TF=="IFITM3", target %in% plot.genes) %>%
                   add_row(TF="IFITM3", target="IFITM3", importance=0),
                 filter(refined.modules, TF=="A2M", target %in% plot.genes) %>%
                   add_row(TF="A2M", target="A2M", importance=0)) %>%
  mutate(y_lab=factor(target, levels=rev(plot.genes)))

p1.1 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8),
			panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))


# plot FTL
plot.genes = filter(refined.modules, TF=="FTL") %>% arrange(desc(importance)) %>% pull(target)
plot.genes = c("FTL", plot.genes)

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))

tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))

tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

p2 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
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

tmp2 = filter(refined.modules, TF=="FTL", target %in% plot.genes) %>% 
	add_row(TF="FTL", target="FTL", importance=0) %>%
  mutate(y_lab=factor(target, levels=rev(plot.genes)))

p2.1 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8),
                        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))


# plot CD74
plot.genes = filter(refined.modules, TF=="CD74") %>% arrange(desc(importance)) %>% pull(target)
plot.genes = c("CD74", plot.genes)

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

p3.1 = ggplot(tmp2, aes(y=y_lab, x=importance))+
  geom_bar(stat="identity", fill="grey")+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8),
                        panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
                        plot.margin = margin(.5,.5,2,0, unit="cm"))

ggsave(file="plots/publication/Figure3/dotplot_DEG-GRNs.pdf",
	marrangeGrob(grobs=list(p1, p1.1, p2, p2.1, p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,2), ncol=7)),
	width=6, height=7)


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
