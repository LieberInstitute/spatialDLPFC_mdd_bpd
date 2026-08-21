setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})
set.seed(123)
cpList <- readRDS("plots/colorPalettes.rds")

# load model results
riskList <- readRDS("processed-data/tmp_PRS_DE/glm-binomial_MDD-Bipolar-PRS_predicting-dxMDD-dxBPD_without-RIN.rda")

# load PRS values
cdata = read.csv("raw-data/PRS/PRS_chosen-p-cutoffs.csv")

colnames(cdata)[8:10] = paste0("prs", substr(colnames(cdata)[8:10], start=0, stop=3))
for(i in colnames(cdata)[8:10]) cdata[,i] = scale(cdata[,i])

cdata$condition = factor(cdata$condition, levels=c("NTC","MDD","BPD"))

# boxplots
p1 <- ggplot(cdata, aes(x=condition, y=prsMDD))+
  ggbeeswarm::geom_beeswarm(aes(color=condition), cex=3, size=.5)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+#scale_shape_manual(values=c("F"=16, "M"=17))+
  geom_boxplot(outliers=F, fill="transparent")+ylim(-3,3)+
  labs(x="diagnosis", y="PRS (z-score)", title="GWAS-MDD")+
  theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), text=element_text(size=6))

p2 <- ggplot(cdata, aes(x=condition, y=prsBPD))+
  ggbeeswarm::geom_beeswarm(aes(color=condition), cex=3, size=.5)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+#scale_shape_manual(values=c("F"=16, "M"=17))+
  geom_boxplot(outliers=F, fill="transparent")+ylim(-3,3)+
  labs(x="diagnosis", y="PRS (z-score)", title="GWAS-BPD")+
  theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), text=element_text(size=6))


# make PRS DFrame long for plotting with facets
cdata2 = bind_rows(filter(cdata, condition %in% c("NTC","MDD")) %>% mutate(predictor="dxMDD"),
                   filter(cdata, condition %in% c("NTC","BPD")) %>% mutate(predictor="dxBPD"),
                   mutate(cdata, predictor="dxAny")) %>%
  mutate(predictor=factor(predictor, levels=c("dxMDD","dxBPD","dxAny"),
                          labels=c("NTC+MDD\n(n=79)","NTC+BD\n(n=80)","NTC+DX (MDD or BD)\n(n=119)")))


plist <- lapply(names(riskList), function(x) {
  # extract model p value
  rlist = riskList[[x]]
  pvals = sapply(rlist, function(y) {
    signif(coef(summary(y))["PRS_scaled", "Pr(>|z|)"], 2)
  })
  p.df = data.frame("predictor"=factor(names(pvals), levels=c("dxMDD","dxBPD","dxAny"),
                                       labels=c("NTC+MDD\n(n=79)","NTC+BD\n(n=80)","NTC+DX (MDD or BD)\n(n=119)")), 
                    "pval"=paste("p =",pvals))
  
  # PRS
  cdata2$plot_me = cdata2[[x]]
  
  # fitted values
  tmp = do.call(rbind, lapply(rlist, function(y) {
    tmp = data.frame("brnum"=y$data$brnum, "condition2"=y$data$condition, "fitted"=y$fitted.values)
    tmp$condition2 = droplevels(factor(tmp$condition2, levels=c("NTC","MDD","BPD","DX")))
    tmp$predictor = factor(paste(levels(tmp$condition2), collapse=" "), levels=c("NTC MDD","NTC BPD","NTC DX"),
                           labels=c("NTC+MDD\n(n=79)","NTC+BD\n(n=80)","NTC+DX (MDD or BD)\n(n=119)"))
    return(tmp)
  })) 

  df1 = left_join(cdata2, tmp)

  ggplot(df1, aes(x=plot_me, y=fitted))+
    geom_point(aes(color=condition2), size=.3)+scale_color_manual(values=c(cpList$dx.pal,"DX"="black"))+
    geom_text(data=p.df, aes(x=-3, y=0, label=pval), hjust=0, vjust=0, size=2)+
    geom_hline(aes(yintercept=.5), lty=2)+ylim(0,1)+xlim(-3,3)+
    facet_wrap(vars(predictor))+
    labs(x=paste0(x, " (z-score)"), y="predicted diagnosis")+
    theme_bw()+theme(text=element_text(size=6), legend.key.size=unit(6,"pt"),
                     panel.grid.minor=element_blank())
})

ggsave(file="plots/publication/Figure_eQTL/supp_dx-PRS-model.pdf",
       grid.arrange(p1, p2, plist[[1]], plist[[2]], layout_matrix=rbind(c(1,3,3,3), c(2,4,4,4))),
       width=6, height=3)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
