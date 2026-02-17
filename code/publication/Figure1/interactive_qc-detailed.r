library(dplyr)
library(ggplot2)

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
colnames(cdata)
table(cdata[,c("remove_problem.areas","lowumi","in_tissue")])

#remove off tissue spots and extra MDD M sample
cdata = filter(cdata, in_tissue==T, sum_umi>0, sample_id!="V13F27-338_C1") %>%
  mutate(qc_filter = ifelse(remove_spots==F, "keep", NA))

cdata[cdata$spotsweeper_outlier, "qc_filter"] = "SpotSweeper"
cdata[cdata$true_edges, "qc_filter"] = "edge artifact"
cdata[cdata$lowumi & is.na(cdata$qc_filter), "qc_filter"] = "<100 UMI"
cdata[cdata$problem_areas_binary & is.na(cdata$qc_filter),"qc_filter"] = "problem area"
table(cdata$qc_filter, useNA="ifany")
cdata$qc_filter = factor(cdata$qc_filter, levels=c("keep","SpotSweeper","edge artifact","<100 UMI","problem area"))

cdata$qc_pass = factor(cdata$qc_filter=="keep", levels=c(TRUE, FALSE), labels=c("keep","discard"))

cdata$condition = factor(cdata$condition, levels=c("NTC","MDD","BPD"))


cdata.n = group_by(cdata, condition, sex) %>% add_tally(name="total_spots") %>%
  group_by(condition, sex, total_spots, qc_filter) %>% tally(name="n_spots") %>%
  mutate(prop_spots=round((n_spots/total_spots)*100, 2)) %>% 
  select(condition, sex, total_spots, qc_filter, prop_spots) %>%
  tidyr::pivot_wider(names_from="qc_filter", values_from="prop_spots")

gt1 <- tableGrob(cdata.n[,4:8], rows=paste(cdata.n$condition, cdata.n$sex))




p1 <- ggplot(cdata, aes(x=condition, y=sum_umi, fill=qc_pass))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(y="sum_umi (log10 scale)", fill="", title="Library size")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

p2 <- ggplot(cdata, aes(x=condition, y=sum_gene, fill=qc_pass))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:4)), labels=c("1","10","100","1k","10k"))+
  labs(y="sum_gene (log10 scale)", fill="", title="Detected genes")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

p3 <- ggplot(cdata, aes(x=condition, y=expr_chrM_ratio, fill=qc_pass))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  labs(y="expr_chrM_ratio", fill="", title="Mitochondrial fraction")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

tmp = group_by(cdata, sample_id, condition, sex) %>% add_tally(name="total_spots") %>%
  group_by(sample_id, condition, sex, total_spots, qc_pass) %>% tally() %>%
  mutate(cond_sex=factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")),
         prop_spots=n/total_spots)

p4 <- ggplot(tmp, aes(x=condition, y=prop_spots, color=qc_pass))+
  ggbeeswarm::geom_quasirandom()+
  scale_color_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  labs(y="prop. of in tissue spots", color="", title="Spots passed QC")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))


# split more 

p5 <- ggplot(cdata, aes(x=condition, y=sum_umi, fill=qc_filter))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("white","orange","dodgerblue","tomato","#FFC0B5"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(y="sum_umi (log10 scale)", fill="", title="Library size")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

p6 <- ggplot(cdata, aes(x=condition, y=sum_gene, fill=qc_filter))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("white","orange","dodgerblue","tomato","#FFC0B5"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:4)), labels=c("1","10","100","1k","10k"))+
  labs(x="condition", y="sum_gene (log10 scale)", fill="", title="Detected genes")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

p7 <- ggplot(cdata, aes(x=condition, y=expr_chrM_ratio, fill=qc_filter))+
  geom_boxplot(position=position_dodge2(width=.6), width=.8, outlier.size=.5)+
  scale_fill_manual(values=c("white","orange","dodgerblue","tomato","#FFC0B5"))+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  labs(x="condition", y="expr_chrM_ratio", fill="", title="Mitochondrial fraction")+
  theme_bw()+theme(axis.title.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))

tmp = group_by(cdata, sample_id, condition, sex, qc_filter) %>% tally() %>%
  mutate(cond_sex=factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")))

p8 <- ggplot(tmp, aes(x=sample_id, y=n, fill=qc_filter))+
  geom_bar(stat="identity", position="fill", color="black", linewidth=.3, width=.7)+
  scale_fill_manual(values=c("white","orange","dodgerblue","tomato","#FFC0B5"))+
  facet_wrap(vars(cond_sex), ncol=1, scales="free_x")+
  scale_y_continuous(breaks=c(.25,.5,.75))+
  labs(x="sample_id", y="prop. of in tissue spots", fill="", title="Spots passed QC")+
  theme_bw()+theme(axis.text.x=element_blank(), panel.ontop = T,
                   panel.background=element_blank(), panel.grid.minor=element_blank(),
                   panel.grid.major.x=element_blank(),
                   legend.box.spacing = unit(0,"pt"), 
                   legend.background=element_rect(fill="transparent",color="transparent"))


pdf(file="plots/publication/Figure1/qc_detailed.pdf")
grid.arrange(p1, p2, p3, p4, ncol=2)
grid.arrange(gt1, top="Percent of spots kept and discarded, by QC approach:")
p5
p6
p7
p8
dev.off()


