library(SpatialExperiment)
library(HDF5Array)
library(dplyr)
library(ggplot2)
library(escheR)
library(ggrastr)
library(gridExtra)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("\nDim n=120 spe:", dim(spe))

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
cat("\nDim n=120 QC csv:", dim(cdata))


#remove off tissue spots and extra MDD M sample
cdata = filter(cdata, in_tissue==T, sample_id!="V13F27-338_C1") %>%
  mutate(qc_filter = ifelse(remove_spots==F, "keep", NA))
cat("\nDim n=119 QC csv (in tissue spots only):", dim(cdata), "\n")
#548579     41

cdata[cdata$true_edges, "qc_filter"] = "edge artifact"
cdata[cdata$problem_areas_binary & is.na(cdata$qc_filter),"qc_filter"] = "extreme region"
cdata[cdata$lowumi & is.na(cdata$qc_filter), "qc_filter"] = "<100 UMI"
cdata[cdata$spotsweeper_outlier, "qc_filter"] = "SpotSweeper"
cdata$qc_filter = factor(cdata$qc_filter, levels=c("edge artifact","extreme region","<100 UMI","SpotSweeper","keep"),
                         labels=c("extreme","extreme","<=100 UMI","local","keep"))
cat("\nNumber of spots removed and reason:\n")
table(cdata$qc_filter, useNA="ifany")
cat("\nProp. of in tissue spots removed:", round(535248/548579, 4), "\n")
# 0.975699

#subset spe to cdata and add qc filter
spe = spe[,rownames(cdata)]
spe$qc_filter = cdata$qc_filter

stopifnot(identical(colnames(spe), rownames(cdata)))
spe$lowUMI = cdata$lowUMI
spe$true_edges = cdata$true_edges
spe$problem_areas_binary = cdata$problem_areas_binary

spe_sub = spe[,spe$in_tissue & spe$sample_id %in% c("V13Y10-022_C1","V13B23-328_A1","V13B23-301_C1")]
spe_sub$lg10.umi = log10(spe_sub$sum_umi+1)
spe_sub$mt.fract = ifelse(is.na(spe_sub$expr_chrM_ratio), 0, spe_sub$expr_chrM_ratio)

#max(spe_sub$lg10.umi)

tmp = spe_sub[,spe_sub$sample_id=="V13Y10-022_C1"]
tmp = spe_sub[,spe_sub$sample_id=="V13B23-328_A1"]
tmp = spe_sub[,spe_sub$sample_id=="V13B23-301_C1"]

plist <- do.call(c, lapply(c("V13Y10-022_C1","V13B23-328_A1","V13B23-301_C1"), function(x) {
  tmp = spe_sub[,spe_sub$sample_id==x]
  p1 <- plotVisium(tmp, spots=TRUE, image=TRUE, facets=NULL, annotate="lg10.umi", point_size = 1.5)+
    scale_fill_viridis_c(limits=c(0,4.5))
  
  p2 <- plotVisium(tmp, spots=TRUE, image=TRUE, facets=NULL, annotate="qc_filter", point_size = 1)+
    scale_fill_manual(values=c("edge"="tomato","extreme"="tomato","<=100 UMI"="red3","local"="black","keep"="grey"))
  
  list(rasterize(p1, "points", dpi=250), rasterize(p2, "points", dpi=250))
}))

plist2 = lapply(plist, function(x) x+theme(legend.position="none"))
plist2[7:8] = plist[1:2]
ggsave(file="plots/publication/Figure1/qc_spot-plots.pdf",
       marrangeGrob(plist2, ncol=1, nrow=1, top=NULL),
       height=5, width=5)

cat("\nQC spot plots saved to: plots/publication/Figure1/qc_spot-plots.pdf\n")


# boxplots for summary
cdata$condition = factor(cdata$condition, levels=c("NTC","MDD","BPD"))
cdata$qc_pass = factor(cdata$qc_filter=="keep", levels=c(TRUE, FALSE), labels=c("keep","discard"))

p1 <- ggplot(cdata, aes(x=condition, y=sum_umi, fill=qc_pass))+
  geom_boxplot(position=position_dodge(width = .5), width=.8, outliers=F)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(y="sum_umi (log10 scale)", fill="", title="Library size")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(), text=element_text(size=10))


p2 <- ggplot(cdata, aes(x=condition, y=sum_gene, fill=qc_pass))+
  geom_boxplot(position=position_dodge(width = .5), width=.8, outliers=F)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:4)), labels=c("1","10","100","1k","10k"))+
  labs(y="sum_gene (log10 scale)", fill="", title="Detected genes")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(), text=element_text(size=10))

p3 <- ggplot(mutate(cdata, expr_chrM_ratio=ifelse(is.na(expr_chrM_ratio), 0, expr_chrM_ratio)), 
             aes(x=condition, y=expr_chrM_ratio, fill=qc_pass))+
  geom_boxplot(width=.8, outliers=F)+
  coord_cartesian(ylim=c(0,.5))+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  labs(y="expr_chrM_ratio", fill="", title="Mitochondrial fraction")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(), text=element_text(size=10))

tmp = group_by(cdata, sample_id, condition, sex) %>% add_tally(name="total_spots") %>%
  group_by(sample_id, condition, sex, total_spots, qc_pass) %>% tally() %>%
  mutate(prop_spots=n/total_spots)

p4 <- ggplot(tmp, aes(x=condition,  y=prop_spots, color=qc_pass))+
  ggbeeswarm::geom_quasirandom()+
  scale_color_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  labs(y="prop. of in tissue spots", color="", title="Spots passed QC")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(), text=element_text(size=10))




cdata2 = filter(cdata, qc_filter=="keep") %>% 
  mutate(cond_sex= factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")))
order1 = group_by(cdata2, sample_id) %>% summarise(med_gene = median(sum_gene)) %>% 
  arrange(desc(med_gene))

cpList <- readRDS("plots/colorPalettes.rds")
p5 <- ggplot(mutate(cdata2, x_lab= factor(sample_id, levels=order1$sample_id)), 
             aes(x=x_lab, y=sum_gene, fill=condition))+
  geom_boxplot(outlier.size = .1, width=.6)+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(10,10000),
                     breaks=c(10^(1:4)), labels=c("10","100","1k","10k"))+
  scale_fill_manual(values=cpList$dx.pal)+
  facet_wrap(vars(cond_sex), ncol=2, scales="free_x")+
  labs(subtitle="Spots that passed QC: Detected genes", x="sample ID", y="sum_gene (log10 scale)")+
  theme_minimal()+theme(axis.text.x= element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank())

p6 <- ggplot(mutate(cdata2, x_lab= factor(sample_id, levels=order1$sample_id)), 
             aes(x=x_lab, y=sum_gene, fill=condition))+
  geom_violin(scale = "width", color="transparent")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(10,10000),
                     breaks=c(10^(1:4)), labels=c("10","100","1k","10k"))+
  scale_fill_manual(values=cpList$dx.pal)+
  facet_wrap(vars(cond_sex), ncol=2, scales="free_x")+
  labs(subtitle="Spots that passed QC: Detected genes", x="sample ID", y="sum_gene (log10 scale)")+
  theme_minimal()+theme(axis.text.x= element_blank(), panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(), panel.ontop = T
                        )

pdf(file="plots/publication/Figure1/qc_summary-boxplots.pdf",width=6, height=9)
grid.arrange(p1, p2, p3, p4, ncol=1)
p5
p6
dev.off()

cat("\nQC spot plots saved to: plots/publication/Figure1/qc_summary-boxplots.pdf\n")


