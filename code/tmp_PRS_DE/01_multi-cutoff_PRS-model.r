setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)
cpList <- readRDS("plots/colorPalettes.rds")

source("code/tmp_PRS_DE/01_custom-functions.r")

# load for plottings
mdd.df = prs.df = read.csv(paste0("raw-data/PRS/PRS_","MDD",".csv"))
bd.df = prs.df = read.csv(paste0("raw-data/PRS/PRS_","Bipolar",".csv"))
stopifnot(identical(mdd.df$IID, bd.df$IID))
cdata = read.csv("raw-data/sample_info/DLPFC_cross-disorders_demographics_MBv.csv")

cat("\nTop thresholded PRS correlations...\n")
round(cor(mdd.df[,2:5], bd.df[,2:5]), 3)

# logistic regression model for prsMDD
prsMDD = prs_wrapper("MDD")
saveRDS(prsMDD, "processed-data/tmp_PRS_DE/glm-binomial_MDD-PRS_predicting-dxMDD-dxBPD-dxAny.rda")
pmtx = sapply(prsMDD, getPval)
cat("\nprsMDD...\n")
print(pmtx)
apply(pmtx, MARGIN=1, FUN=p.adjust, method="BH")

# plot prsMDD
tmp = cbind("IID"=mdd.df[,1], as.data.frame(scale(mdd.df[,2:5])))
prs.df = left_join(tmp, cdata, by=c("IID"="brnum")) %>%
  mutate(sex=factor(sex, levels=c("F","M")),
         condition=factor(condition, levels=c("NTC","MDD","BPD"))) %>%
  tidyr::pivot_longer(colnames(tmp)[2:5], names_to="cutoff", values_to="PRS_scaled")

p1 <- ggplot(prs.df, aes(x=condition, y=PRS_scaled, fill=condition))+
  geom_boxplot(outlier.size=.5)+facet_grid(cols=vars(cutoff))+
  ylim(-4,4)+scale_fill_manual(values=cpList$dx.pal, guide="none")+
  theme_minimal()+labs(x="", y="PRS (z-score)", title="GWAS-MDD")+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), text=element_text(size=6))


# logistic regression model for prsBPD
prsBPD = prs_wrapper("Bipolar")
saveRDS(prsBPD, "processed-data/tmp_PRS_DE/glm-binomial_BPD-PRS_predicting-dxMDD-dxBPD-dxAny.rda")
pmtx = sapply(prsBPD, getPval)
cat("\nprsBPD...\n")
print(pmtx)
apply(pmtx, MARGIN=1, FUN=p.adjust, method="BH")

# plot prsBPD
tmp = cbind("IID"=bd.df[,1], as.data.frame(scale(bd.df[,2:5])))
prs.df = left_join(tmp, cdata, by=c("IID"="brnum")) %>%
  mutate(sex=factor(sex, levels=c("F","M")),
         condition=factor(condition, levels=c("NTC","MDD","BPD"))) %>%
  tidyr::pivot_longer(colnames(tmp)[2:5], names_to="cutoff", values_to="PRS_scaled")

p2 <- ggplot(prs.df, aes(x=condition, y=PRS_scaled, fill=condition))+
  geom_boxplot(outlier.size=.5)+facet_grid(cols=vars(cutoff))+
  ylim(-4,4)+scale_fill_manual(values=cpList$dx.pal, guide="none")+
  theme_minimal()+labs(x="", y="PRS (z-score)", title="GWAS-BPD")+
  theme(panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(), text=element_text(size=6))

ggsave(file="plots/tmp_PRS_DE/top-thresholds.pdf", grid.arrange(p1, p2, ncol=1),
	width=5, height=3)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

