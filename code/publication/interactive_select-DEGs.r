library(SpatialExperiment)
library(dplyr)
library(ggplot2)

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")

spe_pseudo$seurat_label2 = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Oligo","Astro",
                                                                    "L2.3","L4","L5","L6","Inhb"),
                                  labels=c("M.V","Oligo","Astro",
                                           "L2.3","L4","L5","L6","Inhb"))

adj.results = read.csv("processed-data/07_dx_DE/layer-adjusted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))



plot.genes = c("APOLD1","ELK1")
for (j in plot.genes) {
  colData(spe_pseudo)[[j]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

plist[[4]] <- ggplot(summ.df, aes(x=condition, y=logcounts, color=condition))+
  ggbeeswarm::geom_quasirandom(width=.4)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  scale_color_manual(values=cpList$dx.pal, guide="none")+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  facet_grid(cols=vars(sex), rows=vars(key_genes), scales="free_y", switch="y",
             labeller= as_labeller(c("F"="Female","M"="Male","APOLD1"="APOLD1","ELK1"="ELK1")))+
  labs(x="")+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                 strip.text.x= element_text(size=12), axis.text.x =element_text(size=10),
                 strip.text.y.left = element_text(angle=0, size=12, face="italic"), 
                 aspect.ratio=1)

#ggsave(file="plots/publication/violin.png", p3,
#       bg="white", height=7, width=7)


filter(adj.results, gene_name=="ELK1") %>% select(logFC, adj.P.Val, gene_name, coef, sex)
filter(adj.results, gene_name=="APOLD1") %>% select(logFC, adj.P.Val, gene_name, coef, sex)


plot.genes = c("SST","CORT","CRH")
for (j in plot.genes) {
  colData(spe_pseudo)[[j]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

plist[[5]] <- ggplot(summ.df, aes(x=condition, y=logcounts, color=condition))+
  ggbeeswarm::geom_quasirandom(width=.4)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  scale_color_manual(values=cpList$dx.pal, guide="none")+
  scale_y_continuous(position="right", expand=expansion(mult = c(.05, .1)))+
  facet_grid(cols=vars(sex), rows=vars(key_genes), scales="free_y", switch="y",
             labeller= as_labeller(c("F"="Female","M"="Male","SST"="SST","CORT"="CORT", "CRH"="CRH")))+
  labs(x="")+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                   strip.text.x= element_text(size=12),
                   strip.text.y.left = element_text(angle=0, size=12, face="italic"), 
                   aspect.ratio=1)

ggsave(file="plots/publication/select_DEGs.pdf", 
       gridExtra::marrangeGrob(plist, ncol=1, nrow=1, top=NULL),
       height=7, width=7)

filter(adj.results, gene_name=="SST") %>% select(logFC, adj.P.Val, gene_name, coef, sex)
filter(adj.results, gene_name=="CORT") %>% select(logFC, adj.P.Val, gene_name, coef, sex)
filter(adj.results, gene_name=="CRH") %>% select(logFC, adj.P.Val, gene_name, coef, sex)
