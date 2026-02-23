library(dplyr)
library(ggplot2)
set.seed(123)


cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

results_set="smoothed-k9-1663"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         sex.group=factor(coef, levels=comparisons2),
         cluster="L-A") 

lat.filt = filter(lat, gene_id %in% la.degs$gene_id) %>% mutate(is_DEG=adj.P.Val<.05)
summary(lat.filt$t)

filter(lat.filt, is_DEG) %>% group_by(dir, sex.group) %>% slice_min(n=1, abs(t))
#3.67 is the t cutoff for significance

all.3 = filter(lat.filt, sex.group %in% c("F_NTC.MDD","F_NTC.BPD","M_NTC.BPD")) %>% 
  group_by(gene_id, gene_name) %>% summarise(is_DEG=sum(is_DEG)) %>%
  filter(is_DEG==3) %>% pull(gene_name)
length(all.3) #11

plot.genes = all.3
check = c("CORT","SST","CRH",
          "TNFSF10","ABCG2","A2M","RERGL","RBM3","APOLD1",
          "TEF","ELK1")
lrt = read.csv(pastechecklrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$smoothed.bright)),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2))

tmp = bind_rows(filter(lat, gene_name %in% plot.genes) %>%
                  # change adj.P.Val for genes that don't have sig L-R interaction term
                  mutate(adj.P.Val= ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5)),
                filter(lrt, gene_name %in% plot.genes) %>%
                  # change adj.P.Val for genes that don't have sig L-R interaction term
                  mutate(adj.P.Val= ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5)))

tmp$adj.P.Val_bin = as.numeric(as.character(cut(tmp$adj.P.Val, breaks=c(0,.0001,.01,.05,1), labels=c(4,3,2,1))))
table(tmp$adj.P.Val_bin, useNA="ifany")

summary(tmp$logFC)

tmp$plot.genes = factor(tmp$gene_name, levels=rev(check))
tmp$cluster = factor(tmp$cluster, levels=c("L-A",levels(lrt$cluster)))

p1 <- ggplot(tmp, aes(x=cluster, y=plot.genes, fill=logFC, size=adj.P.Val_bin))+
  geom_count(shape=21, color="grey")+
  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
                       limits=c(-2,2))+
  facet_grid(cols=vars(sex.group))+
  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
  scale_size_continuous("adj. p",limits=c(1,4), range=c(1,6), labels=c(">.05","<.05","<.01","<.0001"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8),
                   legend.position="bottom", axis.title.y=element_blank(),
                   legend.text = element_text(size=8))


mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv")
col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)

tmp2 = filter(mratio.sn, gene_name %in% tmp$gene_name) %>% 
  mutate(gene_name=factor(gene_name, levels=levels(tmp$plot.genes)),
         cellType.target= factor(cellType.target, levels=c("InhN","ExcN","Astro","Oligo","Micro","Vasc")))

p2 = ggplot(tmp2, aes(y=gene_name, x=mean.target, fill=cellType.target))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=col.pal)+
  labs(title=" ")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        legend.position="bottom", legend.title=element_blank(),
                        legend.key.size = unit(10,"pt"))


#gridExtra::grid.arrange(p2, p2.1, ncol=2)
#gridExtra::grid.arrange(p1,p2, layout_matrix=matrix(c(1,1,1,1,1,2), ncol=6))

load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")
spe_summ
gids = unique(tmp$gene_id)
df1 = as.data.frame(assay(spe_summ, "logcounts.prop.detected")[gids,])
tmp3 = tidyr::pivot_longer(tibble::rownames_to_column(df1, var="gene_id"), all_of(colnames(df1)), names_to="sample", values_to="prop.spots.detected") %>%
  left_join(as.data.frame(rowData(spe_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(gene_name= factor(gene_name, levels=levels(tmp$plot.genes)))

p3 <- ggplot(tmp3, aes(y=gene_name, x=prop.spots.detected))+
  geom_boxplot(outlier.size=.5, fill="grey")+xlim(0,1)+
  #scale_fill_manual(values=cpList$smoothed.light)+
  labs(title=" ", x="prop. spots\ndetected")+
  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        plot.margin = margin(.2,0,1.5,0,"cm"))


#gridExtra::grid.arrange(p1,p2, p3, layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7))


ggsave(file="plots/publication/Figure2/dotplot_LA-LR-DEGs_all-3-sig.pdf",
       gridExtra::arrangeGrob(grobs=list(p1,p2, p3), layout_matrix=matrix(c(1,1,1,1,1,2,3), ncol=7)),
       width=6, height=5)
