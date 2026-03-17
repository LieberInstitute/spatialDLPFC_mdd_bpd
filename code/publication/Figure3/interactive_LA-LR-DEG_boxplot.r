library(SpatialExperiment)
library(dplyr)
library(ggplot2)

cpList <- readRDS("plots/colorPalettes.rds")
cpList$transfer.bright = cpList$transfer.bright[c(1:4,8,5:7)]

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")

results_set="smoothed-k9-1663"
la.degs_sm = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
lr.degs_sm = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
sm.degs = union(la.degs_sm$gene_name, lr.degs_sm$gene_name)
length(sm.degs) #649

sm.degs2 = union(la.degs_sm$gene_name[la.degs_sm$n_ttest_sig>0], lr.degs_sm$gene_name[lr.degs_sm$n_ttest_sig>0])
length(sm.degs2) #387

results_set="seurat-pc30"
la.degs_se = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
lr.degs_se = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
se.degs = union(la.degs_se$gene_name, lr.degs_se$gene_name)
length(se.degs) #627

se.degs2 = union(la.degs_se$gene_name[la.degs_se$n_ttest_sig>0], lr.degs_se$gene_name[lr.degs_se$n_ttest_sig>0])
length(se.degs2) #408

mbv.degs = sort(union(sm.degs, se.degs))
length(mbv.degs) #817

mbv.degs2 = sort(union(sm.degs2, se.degs2))
length(mbv.degs2) #503



plot.genes = c("SST","CORT","CRH")

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% 
  mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright))))

names(plot.genes) <- plot.genes

all.df_filt = filter(all.df, key_genes %in% plot.genes) 
all.df_filt2 = group_by(all.df_filt, sex, key_genes) %>% summarise(ypos=max(logcounts)) %>%
  group_by(key_genes) %>% mutate(ypos=max(ypos)+.1) %>% arrange(key_genes)

p2 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
  geom_boxplot(aes(fill=condition), color="black", position = position_dodge2(width=.8), 
               outliers=F, 
               linewidth=.2)+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  scale_x_discrete(labels=c("L-A","M.V","Ast","L2.3","L4","In","L5","L6","Olg"))+
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

p3 <- p2+geom_text(data=all.df_filt2, aes(x="L4", y=ypos, label=key_genes), 
                   color="grey50", size=2, fontface="italic", hjust=.5)+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  theme(strip.text.y.left = element_blank())

ggsave(file="plots/publication/Figure3/deg-boxplot_SST-CORT-CRH.pdf", p3, height=3, width=3)

filter(lr.degs_se, gene_name=="SST")
filter(lr.degs_se, gene_name=="CORT")



filter(la.degs_se, gene_name %in% c("C3","CX3CR1","LAPTM5","CSF1R","HLA-DPA1"))
#M_NTC.BPD sig for all
filter(lr.degs_se, gene_name %in% c("C3","CX3CR1","LAPTM5","CSF1R","HLA-DPA1"))
#L6, L5, L2.3 for CX3CR1 and C3, also Astro for CX3CR1


filter(la.degs_se, gene_name %in% c("C3","CX3CR1","LAPTM5","CSF1R","HLA-DPA1"))
#M_NTC.BPD sig for all
filter(lr.degs_se, gene_name %in% c("GADD45B","IFITM2","IFITM3"))
#L6, L5, L2.3 for CX3CR1 and C3, also Astro for CX3CR1

filter(la.degs_se, gene_name %in% c("C1QB","FKBP5"))
filter(lr.degs_se, gene_name %in% c("S100A10","SELPLG","SOCS3"))
sort(filter(lr.degs_se, n_ttest_sig_M_NTC.BPD>0, med_prop.detected_)$gene_name)
#S100A10

plot.genes = c("CX3CR1","C3","GADD45B","IFITM3")

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% 
  mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright))))

names(plot.genes) <- plot.genes

all.df_filt = filter(all.df, key_genes %in% plot.genes) 
all.df_filt2 = group_by(all.df_filt, sex, key_genes) %>% summarise(ypos=max(logcounts)) %>%
  group_by(key_genes) %>% mutate(ypos=max(ypos)+.1) %>% arrange(key_genes)

p2 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
  geom_boxplot(aes(fill=condition), color="black", position = position_dodge2(width=.8), 
               outliers=F, 
               linewidth=.2)+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  scale_x_discrete(labels=c("L-A","M.V","Ast","L2.3","L4","In","L5","L6","Olg"))+
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

p3 <- p2+geom_text(data=all.df_filt2, aes(x="L4", y=ypos, label=key_genes), 
                   color="grey50", size=2, fontface="italic", hjust=.5)+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  theme(strip.text.y.left = element_blank())

ggsave(file="plots/publication/Figure3/deg-boxplot_CX3CR1-C3-GADD45B-IFITM3.pdf", p3, height=4, width=3)

filter(la.degs_se, gene_name %in% c("GADD45B","IFITM3"))

check = filter(lr.degs_se, n_ttest_sig_Astro>0, med_prop.spots.detected_Astro>.03)$gene_name
mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv")
mratio.sn2 = group_by(mratio.sn, gene) %>% slice_min(n=1, MeanRatio.rank)
filter(mratio.sn2, gene_name %in% check, cellType.target=="Astro", MeanRatio>1.5) %>%
  select(gene_name, cellType.target, mean.target, MeanRatio, prop.detected)
#only 3 genes: VEGFA, ITGB4, ANGPTL4

plot.genes = c("ABCG2","VEGFA","ITGB4","ANGPTL4")

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% 
  mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright))))

names(plot.genes) <- plot.genes

all.df_filt = filter(all.df, key_genes %in% plot.genes) 
all.df_filt2 = group_by(all.df_filt, sex, key_genes) %>% summarise(ypos=max(logcounts)) %>%
  group_by(key_genes) %>% mutate(ypos=max(ypos)+.1) %>% arrange(key_genes)

p2 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
  geom_boxplot(aes(fill=condition), color="black", position = position_dodge2(width=.8), 
               outliers=F, 
               linewidth=.2)+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  scale_x_discrete(labels=c("L-A","M.V","Ast","L2.3","L4","In","L5","L6","Olg"))+
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

p3 <- p2+geom_text(data=all.df_filt2, aes(x="L4", y=ypos, label=key_genes), 
                   color="grey50", size=2, fontface="italic", hjust=.5)+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  theme(strip.text.y.left = element_blank())

ggsave(file="plots/publication/Figure3/deg-boxplot_ABCG2-VEGFA-ITGB4-ANGPTL4.pdf", p3, height=4, width=3)


filter(la.degs_se, gene_name %in% plot.genes)
filter(lr.degs_se, gene_name %in% plot.genes)
filter(lr.degs_se, gene_name %in% c("CRH"))




check = filter(lr.degs_se, n_ttest_sig_Micro.Vasc>0, med_prop.spots.detected_Micro.Vasc>.03)$gene_name

check2 = filter(mratio.sn2, gene_name %in% check, cellType.target=="Vasc", MeanRatio>1.5) %>%
  select(gene_name, cellType.target, mean.target, MeanRatio, prop.detected) %>%
  pull(gene_name)

setdiff(check2, sm.degs)
