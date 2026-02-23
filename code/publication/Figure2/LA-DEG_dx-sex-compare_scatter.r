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

results_set = "smoothed-k9-1663"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
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

tmp1 = filter(lat.filt, sex.group %in% c("F_NTC.MDD","F_NTC.BPD")) %>%
                   select(gene_id, gene_name, sex.group, t) %>% 
                   tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(F_NTC.MDD)>=3.67, abs(F_NTC.BPD)>=3.67),
         dir=paste(sign(F_NTC.MDD), sign(F_NTC.BPD)))
tmp1$pt.col = ifelse(tmp1$sig=="FALSE FALSE", "grey", "black")
tmp1[tmp1$sig=="TRUE TRUE" & tmp1$dir=="1 1","pt.col"] = "red3"
tmp1[tmp1$sig %in% c("TRUE FALSE") & tmp1$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp1[tmp1$sig %in% c("FALSE TRUE") & tmp1$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp1[tmp1$sig=="TRUE TRUE" & tmp1$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp1[tmp1$sig %in% c("TRUE FALSE") & tmp1$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp1[tmp1$sig %in% c("FALSE TRUE") & tmp1$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes1 = c(all.3, c("RAMP2","VGF","SURF1"), c("DDIT4","HBB"),
                c("RASD1","MT2A","HILPDA","CEBPD","GADD45B","ANGPTL4","JUN","FCER1G"),
                c("FOXO3","SPOP","MAPK3"),
                c("CIRBP","LMO2","SCG2","CLDN11","MAG"),
                c("ATP6V0E2"),
                c("MUSTN1"), c("SRGN","FCGR2A"), c("C1QB"),
                c("CNDP1","RNASE1"))

p1 <- ggplot(tmp1, aes(x=F_NTC.MDD, y=F_NTC.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="MDD and BPD common (vs NTC) (Females)", y="F_NTC.BPD t stat", x="F_NTC.MDD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p1t <- p1+ggrepel::geom_text_repel(data=filter(tmp1, gene_name %in% plot.genes1) %>%
                             mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                  labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                           aes(label=gene_name),
                           min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp2 = filter(lat.filt, sex.group %in% c("M_NTC.MDD","M_NTC.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(M_NTC.MDD)>=3.67, abs(M_NTC.BPD)>=3.67),
         dir=paste(sign(M_NTC.MDD), sign(M_NTC.BPD)))
tmp2$pt.col = ifelse(tmp2$sig=="FALSE FALSE", "grey", "black")
tmp2[tmp2$sig=="TRUE TRUE" & tmp2$dir=="1 1","pt.col"] = "red3"
tmp2[tmp2$sig %in% c("TRUE FALSE") & tmp2$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp2[tmp2$sig %in% c("FALSE TRUE") & tmp2$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp2[tmp2$sig=="TRUE TRUE" & tmp2$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp2[tmp2$sig %in% c("TRUE FALSE") & tmp2$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp2[tmp2$sig %in% c("FALSE TRUE") & tmp2$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes2 = c(all.3, 
                c("CIRBP", "MAPK3", "ZFP36", "EDN1", "PDLIM4", "ANGPTL4", "VASN", "BAIAP3", "GADD45B", "CEBPD"),
                c("HLA-DPA1","CSF1R","C3","FOLR2", "CX3CR1", "DLX6-AS1"),
                c("ALDOA","EDNRB"), c("ZNF385D"))

p2 <- ggplot(tmp2, aes(x=M_NTC.MDD, y=M_NTC.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="MDD and BPD common (vs NTC) (Males)", y="M_NTC.BPD t stat", x="M_NTC.MDD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p2t <- p2+ggrepel::geom_text_repel(data=filter(tmp2, gene_name %in% plot.genes2) %>%
                                     mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                          labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                   aes(label=gene_name),
                                   min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp3 = filter(lat.filt, sex.group %in% c("F_NTC.MDD","F_MDD.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(F_NTC.MDD)>=3.67, abs(F_MDD.BPD)>=3.67),
         dir=paste(sign(F_NTC.MDD), sign(-F_MDD.BPD)))
tmp3$pt.col = ifelse(tmp3$sig=="FALSE FALSE", "grey", "black")
tmp3[tmp3$sig=="TRUE TRUE" & tmp3$dir=="1 1","pt.col"] = "red3"
tmp3[tmp3$sig %in% c("TRUE FALSE") & tmp3$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp3[tmp3$sig %in% c("FALSE TRUE") & tmp3$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp3[tmp3$sig=="TRUE TRUE" & tmp3$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp3[tmp3$sig %in% c("TRUE FALSE") & tmp3$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp3[tmp3$sig %in% c("FALSE TRUE") & tmp3$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes3 = c(all.3, c("RASD1","HILPDA","CEBPD","DDIT4","GADD45B","ANGPTL4","SESN1","JUN"),
                c("FCGR2A","FCER1G"),
                c("SURF1"),
                c("SRGN","RGS1"), c("SPP1", "C1QB"),
                c("TF","CNDP1","CLDN11"),
                c("MUSTN1"), c("TTYH1","CRTC1"))

p3 <- ggplot(tmp3, aes(x=F_NTC.MDD, y=F_MDD.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="MDD common (vs NTC and BPD) (Females)", x="F_NTC.MDD t stat", y="F_MDD.BPD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p3t <- p3+ggrepel::geom_text_repel(data=filter(tmp3, gene_name %in% plot.genes3) %>%
                                  mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                       labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                aes(label=gene_name),
                                min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp4 = filter(lat.filt, sex.group %in% c("F_NTC.BPD","F_MDD.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(F_NTC.BPD)>=3.67, abs(F_MDD.BPD)>=3.67),
         dir=paste(sign(F_NTC.BPD), sign(F_MDD.BPD)))
tmp4$pt.col = ifelse(tmp4$sig=="FALSE FALSE", "grey", "black")
tmp4[tmp4$sig=="TRUE TRUE" & tmp4$dir=="1 1","pt.col"] = "red3"
tmp4[tmp4$sig %in% c("TRUE FALSE") & tmp4$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp4[tmp4$sig %in% c("FALSE TRUE") & tmp4$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp4[tmp4$sig=="TRUE TRUE" & tmp4$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp4[tmp4$sig %in% c("TRUE FALSE") & tmp4$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp4[tmp4$sig %in% c("FALSE TRUE") & tmp4$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"


plot.genes4 = c(all.3, c("HBB","FOXO3","SPOP","MAPK3","DDIT4"),
                c("SURF1","LMO2"),
                c("CNDP1","CLDN11","MAG"),
                c("RASD1","SRGN"), c("SPP1","C1QB","FCER1G","FCGR2A"),
                c("TTYH1","CRTC1"))

p4 <- ggplot(tmp4, aes(x=F_NTC.BPD, y=F_MDD.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="BPD common (vs NTC and MDD) (Females)", x="F_NTC.BPD t stat", y="F_MDD.BPD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p4t <- p4+ggrepel::geom_text_repel(data=filter(tmp4, gene_name %in% plot.genes4) %>%
                               mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                    labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                             aes(label=gene_name),
                             min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp8 = filter(lat.filt, sex.group %in% c("M_NTC.MDD","M_MDD.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(M_NTC.MDD)>=3.67, abs(M_MDD.BPD)>=3.67),
         dir=paste(sign(M_NTC.MDD), sign(-M_MDD.BPD)))
tmp8$pt.col = ifelse(tmp8$sig=="FALSE FALSE", "grey", "black")
tmp8[tmp8$sig=="TRUE TRUE" & tmp8$dir=="1 1","pt.col"] = "red3"
tmp8[tmp8$sig %in% c("TRUE FALSE") & tmp8$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp8[tmp8$sig %in% c("FALSE TRUE") & tmp8$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp8[tmp8$sig=="TRUE TRUE" & tmp8$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp8[tmp8$sig %in% c("TRUE FALSE") & tmp8$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp8[tmp8$sig %in% c("FALSE TRUE") & tmp8$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes8 = c(all.3, c("EDNRB","ALDOA","CIRBP"),
                c("DLX6-AS1","GADD45B","PDLIM4","EDN1","APOLD1","ANGPTL4"),
                c("COX7A1","ZFP36","ZFP36L2"))

p8 <- ggplot(tmp8, aes(x=M_NTC.MDD, y=M_MDD.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="MDD common (vs NTC and BPD) (Males)", x="M_NTC.MDD t stat", y="M_MDD.BPD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p8t <- p8+ggrepel::geom_text_repel(data=filter(tmp8, gene_name %in% plot.genes8) %>%
                                     mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                          labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                   aes(label=gene_name),
                                   min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp5 = filter(lat.filt, sex.group %in% c("M_NTC.BPD","M_MDD.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(M_NTC.BPD)>=3.67, abs(M_MDD.BPD)>=3.67),
         dir=paste(sign(M_NTC.BPD), sign(M_MDD.BPD)))
tmp5$pt.col = ifelse(tmp5$sig=="FALSE FALSE", "grey", "black")
tmp5[tmp5$sig=="TRUE TRUE" & tmp5$dir=="1 1","pt.col"] = "red3"
tmp5[tmp5$sig %in% c("TRUE FALSE") & tmp5$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp5[tmp5$sig %in% c("FALSE TRUE") & tmp5$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp5[tmp5$sig=="TRUE TRUE" & tmp5$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp5[tmp5$sig %in% c("TRUE FALSE") & tmp5$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp5[tmp5$sig %in% c("FALSE TRUE") & tmp5$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes5 = c(all.3, c("GADD45B","ANGPTL4","MAPK3","VASN","CEBPD","TIMP1","EDN1"),
                c("PDLIM4","BAIAP3","CIRBP","GADD45G","ZFP36"),
                c("DLX6-AS1","VIP","GAD2"),
                c("CX3CR1","FOLR2", "CSF1R", "C3", "HLA-DPA1","LAPTM5"),
                c("ZFP36L2","C11orf96"),
                c("COX7A1","ATP5ME"))

p5 <- ggplot(tmp5, aes(x=M_NTC.BPD, y=M_MDD.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="BPD common (vs NTC and MDD) (Males)", x="M_NTC.BPD t stat", y="M_MDD.BPD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p5t <- p5+ggrepel::geom_text_repel(data=filter(tmp5, gene_name %in% plot.genes5) %>%
                                     mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                          labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                   aes(label=gene_name),
                                   min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

tmp6 = filter(lat.filt, sex.group %in% c("M_NTC.BPD","F_NTC.MDD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(M_NTC.BPD)>=3.67, abs(F_NTC.MDD)>=3.67),
         dir=paste(sign(M_NTC.BPD), sign(F_NTC.MDD)))
tmp6$pt.col = ifelse(tmp6$sig=="FALSE FALSE", "grey", "black")
tmp6[tmp6$sig=="TRUE TRUE" & tmp6$dir=="1 1","pt.col"] = "red3"
tmp6[tmp6$sig %in% c("TRUE FALSE") & tmp6$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp6[tmp6$sig %in% c("FALSE TRUE") & tmp6$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp6[tmp6$sig=="TRUE TRUE" & tmp6$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp6[tmp6$sig %in% c("TRUE FALSE") & tmp6$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp6[tmp6$sig %in% c("FALSE TRUE") & tmp6$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes6 = c(all.3, c("RASD1","CEBPD","GADD45B","ELK1","ANGPTL4","ZFP36"),
                c("CIRBP","RBM3","FCGR2A","FCGR3A","FCER1G","BAIAP3"),
                c("RAMP2","VGF","ATP6V0E2","SLC38A5"),
                c("VIP","CX3CR1","GAD2","SLC6A1","DLX6-AS1","FOLR2","PLD4"),
                c("TIMP1","EDN1","PDLIM4","MAPK3","VASN","GADD45G"),
                c("DDIT4","SESN1","BAG3","HSPB1","JUN","FABP5","MT2A","HILPDA"),
                c("SURF1","VEGFA","HLA-DPA1","LAPTM5","CSF1R","C3"))

p6 <- ggplot(tmp6, aes(x=F_NTC.MDD, y=M_NTC.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="MDD females & BPD males common (vs NTC)", y="M_NTC.BPD t stat", x="F_NTC.MDD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p6t <- p6+ggrepel::geom_text_repel(data=filter(tmp6, gene_name %in% plot.genes6) %>%
                                     mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                          labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                   aes(label=gene_name),
                                   min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")


tmp7 = filter(lat.filt, sex.group %in% c("M_NTC.BPD","F_NTC.BPD")) %>%
  select(gene_id, gene_name, sex.group, t) %>% 
  tidyr::pivot_wider(names_from="sex.group", values_from="t") %>%
  mutate(sig=paste(abs(M_NTC.BPD)>=3.67, abs(F_NTC.BPD)>=3.67),
         dir=paste(sign(M_NTC.BPD), sign(F_NTC.BPD)))
tmp7$pt.col = ifelse(tmp7$sig=="FALSE FALSE", "grey", "black")
tmp7[tmp7$sig=="TRUE TRUE" & tmp7$dir=="1 1","pt.col"] = "red3"
#tmp7[tmp7$sig %in% c("TRUE FALSE", "FALSE TRUE") & tmp7$dir=="1 1","pt.col"] = "#FF917E"
tmp7[tmp7$sig %in% c("TRUE FALSE") & tmp7$dir %in% c("1 -1","1 1"),"pt.col"] = "#FF917E"
tmp7[tmp7$sig %in% c("FALSE TRUE") & tmp7$dir %in% c("-1 1","1 1"),"pt.col"] = "#FF917E"

tmp7[tmp7$sig=="TRUE TRUE" & tmp7$dir=="-1 -1","pt.col"] = "dodgerblue3"
tmp7[tmp7$sig %in% c("TRUE FALSE") & tmp7$dir %in% c("-1 -1","-1 1"),"pt.col"] = "skyblue"
tmp7[tmp7$sig %in% c("FALSE TRUE") & tmp7$dir %in% c("-1 -1","1 -1"),"pt.col"] = "skyblue"

plot.genes7 = c(all.3, c("MAPK3","DLX6-AS1"),
                c("CEBPD","GADD45B","ANGPTL4","VASN","GADD45G","CIRBP","PDLIM4"),
                c("SURF1","RAMP2","VGF","SLC38A5"), c("DDIT4"),
                c("CX3CR1","FOLR2","PLD4","HLA-DPA1","CSF1R","C3"),
                c("GAD1","MAG","SCG2","LMO2","CLDN11"),
                c("FOXO3","HBB","FTL","SPOP","PUM1","RPL28"))

p7 <- ggplot(tmp7, aes(x=F_NTC.BPD, y=M_NTC.BPD, color=pt.col))+
  geom_point(size=.5)+scale_color_identity()+
  geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+
  geom_hline(aes(yintercept=-3.67), lty=2)+geom_hline(aes(yintercept=3.67), lty=2)+
  geom_vline(aes(xintercept=-3.67), lty=2)+geom_vline(aes(xintercept=3.67), lty=2)+
  coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
  labs(title="BPD females & BPD males common (vs NTC)", x="F_NTC.BPD t stat", y="M_NTC.BPD t stat")+
  theme_bw()+theme(aspect.ratio=1, text=element_text(size=10), 
                   plot.title.position = "plot", plot.title=element_text(size=9))

p7t <- p7+ggrepel::geom_text_repel(data=filter(tmp7, gene_name %in% plot.genes7) %>%
                                     mutate(pt.col=factor(pt.col, levels=c("dodgerblue3","skyblue","black","#FF917E","red3"),
                                                          labels=c("dodgerblue3","dodgerblue3","black","red3","red3"))), 
                                   aes(label=gene_name),
                                   min.segment.length = 0, max.overlaps = Inf, size=2, fontface="italic")

pdf(file="plots/publication/Figure2/LA-F-test-padj05_dx-sex-compare_scatter-each.pdf", height=3.5, width=3.5)
p1t
p3t
p4t
p2t
p8t
p5t
p7t
p6t
dev.off()


ggsave(file="plots/publication/Figure2/LA-F-test-padj05_dx-sex-compare_scatter-no-genes.pdf", 
       arrangeGrob(grobs=list(p1, p3, p4, 
                              p2, p8, p5,
                              p7, p6), ncol=3), height=9, width=9)
ggsave(file="plots/publication/Figure2/LA-F-test-padj05_dx-sex-compare_scatter-with-genes.pdf", 
       arrangeGrob(grobs=list(p1t, p3t, p4t, 
                              p2t, p8t, p5t,
                              p7t, p6t), ncol=3), height=9, width=9)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()


