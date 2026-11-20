args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript 04_QC-plots.R <qc_df> <color_palettes>")
}

suppressPackageStartupMessages({
    library("dplyr")
    library("ggplot2")
    library("gridExtra")
    library("ggbeeswarm")
})

set.seed(123) # for consistency with ggbeeswarm
cpList <- readRDS(args[[2L]])
qc.df <- read.csv(args[[1L]], row.names=1)

## Format QC annotations for plotting
qc.df <- filter(qc.df, in_tissue)
qc.df$dx = factor(qc.df$dx, levels=c("NTC","MDD","BD"))
qc.df$qc_pass = factor(qc.df$remove_spots, levels=c(FALSE, TRUE), labels=c("keep","discard"))
qc.df$dx_sex= interaction(qc.df$dx, qc.df$sex, sep=" ", drop=TRUE, lex.order=TRUE)



# QC DX-SEX GROUP SUMMARY PLOTS #

# Guard against an implicit default graphics device during plot construction.
pdf(NULL, width=6, height=9)

p1 <- ggplot(qc.df, aes(x=dx, y=sum_umi, fill=qc_pass))+
  geom_boxplot(position=position_dodge(width = .5), width=.8, outliers=F)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
    breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(y="sum_umi (log10 scale)", fill="", title="Library size")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(), text=element_text(size=10))

p2 <- ggplot(qc.df, aes(x=dx, y=sum_gene, fill=qc_pass))+
  geom_boxplot(position=position_dodge(width = .5), width=.8, outliers=F)+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
    breaks=c(10^(0:4)), labels=c("1","10","100","1k","10k"))+
  labs(y="sum_gene (log10 scale)", fill="", title="Detected genes")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(), text=element_text(size=10))

p3 <- ggplot(mutate(qc.df, expr_chrM_ratio=ifelse(is.na(expr_chrM_ratio), 0, expr_chrM_ratio)), 
    aes(x=dx, y=expr_chrM_ratio, fill=qc_pass))+
  geom_boxplot(width=.8, outliers=F)+
  coord_cartesian(ylim=c(0,.5))+
  scale_fill_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  labs(y="expr_chrM_ratio", fill="", title="Mitochondrial fraction")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(), text=element_text(size=10))

tmp = group_by(qc.df, sample_id, dx, sex) %>% add_tally(name="total_spots") %>%
  group_by(sample_id, dx, sex, total_spots, qc_pass) %>% tally() %>%
  mutate(prop_spots=n/total_spots)

p4 <- ggplot(tmp, aes(x=dx,  y=prop_spots, color=qc_pass))+
  ggbeeswarm::geom_quasirandom()+
  scale_color_manual(values=c("grey","tomato"))+
  facet_wrap(vars(sex), ncol=2, scales="free_x")+
  labs(y="prop. of in tissue spots", color="", title="Spots passed QC")+
  theme_minimal()+theme(axis.title.x=element_blank(), panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(), text=element_text(size=10))



# VIOLIN OF DETECTED GENE DISTRIBUTIONS FOR EVERY SAMPLE #

qc.df = filter(qc.df, qc_pass=="keep")
order1 = group_by(qc.df, sample_id) %>% summarise(med_gene = median(sum_gene)) %>% 
  arrange(desc(med_gene))

p5 <- ggplot(mutate(qc.df, x_lab= factor(sample_id, levels=order1$sample_id)), 
    aes(x=x_lab, y=sum_gene, fill=dx))+
  geom_violin(scale = "width", color="transparent")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(10,10000),
    breaks=c(10^(1:4)), labels=c("10","100","1k","10k"))+
  scale_fill_manual(values=cpList$dx.pal)+
  facet_wrap(vars(dx_sex), ncol=2, scales="free_x")+
  labs(subtitle="Spots that passed QC: Detected genes", x="sample ID", y="sum_gene (log10 scale)")+
  theme_minimal()+theme(axis.text.x= element_blank(), panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(), panel.ontop = T
    )


# SAVE QC PLOTS #

dev.off() # end fileless device
pdf(file="qc_summary_plots.pdf",width=6, height=9)
grid.arrange(p1, p2, p3, p4, ncol=1)
print(p5)
dev.off()
