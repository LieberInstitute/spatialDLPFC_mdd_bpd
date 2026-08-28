setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)
cpList <- readRDS("plots/colorPalettes.rds")

# load for plottings
mdd.df = prs.df = read.csv(paste0("raw-data/PRS/PRS_","MDD",".csv"))
bd.df = prs.df = read.csv(paste0("raw-data/PRS/PRS_","Bipolar",".csv"))
stopifnot(identical(mdd.df$IID, bd.df$IID))
cdata = read.csv("raw-data/sample_info/DLPFC_cross-disorders_demographics_MBv.csv")

# plot prsMDD
tmp = cbind("IID"=mdd.df[,1], as.data.frame(scale(mdd.df[,2:5])))
prs.df = left_join(tmp, cdata, by=c("IID"="brnum")) %>%
  mutate(sex=factor(sex, levels=c("F","M")),
         condition=factor(condition, levels=c("NTC","MDD","BPD"))) %>%
  tidyr::pivot_longer(colnames(tmp)[2:5], names_to="cutoff", values_to="PRS_scaled")


p1 <- ggplot(prs.df, aes(x=condition, y=PRS_scaled))+
    geom_boxplot(outliers=F)+
  ggbeeswarm::geom_quasirandom(aes(shape=sex, color=condition), size=.5)+
  facet_wrap(vars(cutoff), ncol=2)+
  ylim(-4,4)+scale_color_manual(values=cpList$dx.pal, guide="none")+
  scale_shape_manual(values=c(19,1))+
  theme_minimal()+labs(x="", y="PRS (z-score)", title="GWAS-MDD")+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), text=element_text(size=6))


ggsave(file="plots/publication/supp_PRS/prsMDD_boxplots.pdf", p1,
	width=3, height=2.5)


# plot prsBPD
tmp = cbind("IID"=bd.df[,1], as.data.frame(scale(bd.df[,2:5])))
prs.df = left_join(tmp, cdata, by=c("IID"="brnum")) %>%
  mutate(sex=factor(sex, levels=c("F","M")),
         condition=factor(condition, levels=c("NTC","MDD","BPD"))) %>%
  tidyr::pivot_longer(colnames(tmp)[2:5], names_to="cutoff", values_to="PRS_scaled")

p2 <- ggplot(prs.df, aes(x=condition, y=PRS_scaled))+
  geom_boxplot(outliers=F)+
  ggbeeswarm::geom_quasirandom(aes(shape=sex, color=condition), size=.5)+
  facet_wrap(vars(cutoff), ncol=2)+
  ylim(-4,4)+scale_color_manual(values=cpList$dx.pal, guide="none")+
  scale_shape_manual(values=c(19,1))+
  theme_minimal()+labs(x="", y="PRS (z-score)", title="GWAS-BPD")+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), text=element_text(size=6))

ggsave(file="plots/publication/supp_PRS/prsBPD_boxplots.pdf", p2,
        width=3, height=2.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
