library(dplyr)
library(ggplot2)
set.seed(123)


cpList <- readRDS("plots/colorPalettes.rds")

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

lat.filt = filter(lat, gene_id %in% la.degs$gene_id) %>% mutate(is_DEG=adj.P.Val<.05)
summary(lat.filt$t)

filter(lat.filt, is_DEG) %>% group_by(dir, sex.group) %>% slice_min(n=1, abs(t))
#3.58 is the t cutoff for significance

all.3 = filter(lat.filt, sex.group %in% c("F_NTC.MDD","F_NTC.BPD","M_NTC.BPD")) %>% 
  group_by(gene_id, gene_name) %>% summarise(is_DEG=sum(is_DEG)) %>%
  filter(is_DEG==3) %>% pull(gene_name)
length(all.3) #16

#check for genes that weren't included in smoothed DE input
sm.la = read.csv("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_rev-gene-input_F-test.csv")
setdiff(all.3, sm.la$gene_name)
#just MT1A, exclude this from main figure


check = c("CORT","SST","CRH",
          "TNFSF10","ABCG2","A2M","RERGL","RBM3","APOLD1",
          "TEF","ELK1")

setdiff(check, all.3) #none
setdiff(all.3, check)
#TRB2, ITM2A, GADD45B, MT1X, MT1A


plot.genes = c("CORT","SST","CRH",
               "TNFSF10","ABCG2","A2M","RERGL","RBM3","TRIB2","ITM2A",
               "APOLD1","GADD45B","MT1X",#"MT1A",
               "TEF","ELK1")

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$transfer.bright)[c(1:4,8,5:7)]),
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

tmp$plot.genes = factor(tmp$gene_name, levels=rev(plot.genes))
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
p1
