library(SpatialExperiment)
library(dplyr)
library(ggplot2)

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")
cpList$transfer.bright = cpList$transfer.bright[c(1:4,8,5:7)]

results_set="smoothed-k9-1663"
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

results_set="seurat-pc30"
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")

la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

plot.genes= c("BCL6","HSPA1B","JUN",#WM-inflam,
              "APLP1","SLC44A1","MAG", #WM-myelin
              "UBA52","SURF1","DDIT4",#trans/mito
              "CEBPD","GADD45B","ANGPTL4", #BBB angio inflamm group
              "SST", "CORT", "CRH", #InhN that are L-R too
              "ELK1","DUSP6","RASD1") #MAPK

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}


summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         smoothed_k9_1663="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "smoothed_k9_1663", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(smoothed_k9_1663, levels=c("L-A", names(cpList$smoothed.bright))))

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright)),
                                                                     labels=c("L-A","M.V","A","L2.3","L4","In","L5","L6","Olg")))


plot.genes= c("BCL6","HSPA1B","JUN") #WM-inflamm,
plot.genes= c("APLP1","SLC44A1","MAG") #WM-myelin
plot.genes= c("UBA52","DDIT4","SURF1") #transl-starv
plot.genes= c("CEBPD","GADD45B","ANGPTL4") #BBB angio inflamm group
plot.genes= c("SST", "CORT", "CRH") #InhN that are L-R too
plot.genes= c("ELK1","DUSP6","RASD1") #MAPK

names(plot.genes) <- plot.genes
all.df_filt = filter(all.df, key_genes %in% plot.genes) 
all.df_filt2 = group_by(all.df_filt, sex, key_genes) %>% summarise(ypos=max(logcounts)) %>%
  group_by(key_genes) %>% mutate(ypos=max(ypos)+.1) %>% arrange(key_genes)

p2 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
  geom_boxplot(aes(fill=condition), color="black", position = position_dodge2(width=.8), 
               outliers=F, 
               linewidth=.2)+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  facet_grid(cols=vars(sex), rows=vars(key_genes), scales="free_y", switch="y",
             labeller= as_labeller(c("F"="Female","M"="Male", plot.genes)))+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                   text=element_text(size=8), axis.text=element_text(size=6),
                   axis.title.x=element_blank(), axis.title.y=element_blank(),
                   axis.ticks = element_line(linewidth=.1),
                   panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.1))

b = ggplot_build(p2)
all.df_filt2 = group_by(b$data[[1]], PANEL) %>% summarise(ym=max(ymax)) %>%
  mutate(key_genes= factor(all.df_filt2$key_genes, levels=levels(all.df_filt2$key_genes))) %>%
  group_by(key_genes) %>% mutate(ypos=max(ym)+.5)

p3 <- p2+geom_text(data=all.df_filt2, aes(x="L3.4", y=ypos, label=key_genes), 
          color="grey50", size=2, fontface="italic", hjust=.5)+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  theme(strip.text.y.left = element_blank())

ggsave(file="plots/publication/Figure2/deg-boxplot_WM-inflamm.pdf", p3, height=3, width=3)
ggsave(file="plots/publication/Figure2/deg-boxplot_WM-myelin.pdf", p3, height=3, width=3)
ggsave(file="plots/publication/Figure2/deg-boxplot_transl-starv.pdf", p3, height=3, width=3)
ggsave(file="plots/publication/Figure2/deg-boxplot_BBB-inflamm.pdf", p3, height=3, width=3)
ggsave(file="plots/publication/Figure2/deg-boxplot_GABA-pep.pdf", p3, height=3, width=3)
ggsave(file="plots/publication/Figure2/deg-boxplot_MAPK.pdf", p3, height=3, width=3)

# for finding sig values to mark with asterisk
filter(la.degs, gene_name %in% plot.genes)
filter(lr.degs, gene_name %in% plot.genes)

# for the dotplot in supp include more genes
plot.genes = c(
               "BCL6","HSPA1B","JUN",#WM-inflam
               "ANP32B","NFKBIA","LMNA", #supp dotplot only
               "APLP1","CYP51A1","KCNMB4","SLC44A1","MAG", #WM-myelin
               #KCNMB4 and myelin: https://pmc.ncbi.nlm.nih.gov/articles/PMC8596180/
               "TF","ENPP2",#supp dotplot only
               "SURF1","ATP6V0E2","DDIT4","SESN1",#trans/mito
               "UBA52","EIF5B","UBC","RPL28","RPS8","RPL29","RPS12",#supp dotplot only
               "CEBPD","GADD45B","ANGPTL4", #BBB angio inflamm group
               "APOLD1","IFITM3","MT2A", #suppdoplot only
               "BAALC-AS1", #supp dotplot only
               "FKBP5", "C1QB","C3","CX3CR1","LAPTM5","CD74","CSF1R","HLA-DPA1", #microglia, supp dotplot only
               "SST", "CORT", "CRH", "VGF","RAMP2",#InhN that are L-R too
               "HSPA8","GPR37", "SLC32A1",#SUPP DOTPLOT ONLY
               "ELK1","MAPK3","DUSP6",
               "DUSP4","RASD1") #supp dotplot only

lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         sex.group=factor(coef, levels=comparisons2),
         cluster="L-A") 

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$smoothed.bright)),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2))

tmp = bind_rows(filter(lat, gene_name %in% plot.genes),
          filter(lrt, gene_name %in% plot.genes))
head(tmp)
tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))
table(tmp$adj.P.Val_bin, useNA="ifany")

summary(tmp$logFC)

tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                                    limits=c(-2,2))+
  facet_grid(cols=vars(sex.group))+
  scale_size_continuous(limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()

