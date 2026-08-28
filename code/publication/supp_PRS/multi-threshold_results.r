setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)
cpList <- readRDS("plots/colorPalettes.rds")

source("code/tmp_PRS_DE/01_custom-functions.r")

# load lists of PRS model results
prsMDD = readRDS("processed-data/tmp_PRS_DE/glm-binomial_MDD-PRS_predicting-dxMDD-dxBPD-dxAny.rda")

prsBPD = readRDS("processed-data/tmp_PRS_DE/glm-binomial_BPD-PRS_predicting-dxMDD-dxBPD-dxAny.rda")

# confidence interval tables
mdd.out = do.call(rbind, lapply(names(prsMDD), function(x) {
  df1 = getConfInt(prsMDD[[x]])
  df1$cutoff = x
  return(df1)
}))

bd.out = do.call(rbind, lapply(names(prsBPD), function(x) {
  df1 = getConfInt(prsBPD[[x]])
  df1$cutoff = x
  return(df1)
}))

all.out = bind_rows(mutate(mdd.out, disorder="MDD"),
          mutate(bd.out, disorder="BPD"))

# p value tables
all.p = bind_rows(mutate(meltPvalDF(prsMDD), disorder="MDD"),
                  mutate(meltPvalDF(prsBPD), disorder="BPD"))

# merge for plotting
all.df = full_join(all.out, all.p, by=c("dx","cutoff","disorder")) %>%
  mutate(disorder=factor(disorder, levels=c("MDD","BPD"), labels=c("prsMDD","prsBPD")),
         dx=factor(dx, levels=rev(c("dxMDD","dxBPD","dxAny"))))

p1 <- ggplot(all.df, aes(x=OR, y=dx))+
  geom_linerange(aes(xmin=CI05, xmax=CI95, group=cutoff, color=adj_p<.05), position=position_dodge(width=.2))+
  scale_color_manual(values=c("FALSE"="grey30", "TRUE"="red"))+
  geom_point(aes(fill=cutoff), position=position_dodge(width=.2), shape=21)+
  scale_fill_manual(values=RColorBrewer::brewer.pal(n=5, "PuBu")[2:5])+
  facet_wrap(vars(disorder), scales="free_y", ncol=1)+
  geom_vline(aes(xintercept=1))+
  theme_bw()+theme(text=element_text(size=6))

ggsave(file="plots/publication/supp_PRS/dx-PRS-model_multi-threshold-results.pdf",
       p1, height=5, width=3)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
