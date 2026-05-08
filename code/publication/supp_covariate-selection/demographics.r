setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")

demo = read.csv("processed-data/publication/demographics.csv")

#chisq
cat("\n\nChi-squared of dx on smoking status...\n")
M = table(demo[,c("condition","Smoking")])
(Xsq <- chisq.test(M))

#lm
cat("\n\nLinear regression of dx*sex on donor covars...")
cat("\n>>> Age...\n")
m1 = lm(age ~ condition*sex, data=demo)
summary(m1)

cat("\n>>> BMI...\n")
m2 = lm(BMI ~ condition*sex, data=demo)
summary(m2)

cat("\n>>> RIN...\n")
m3 = lm(RIN ~ condition*sex, data=demo)
summary(m3)

# now plotting
demo$cond_sex = factor(paste(demo$condition, demo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"),
                       labels=c("NTC\nF","NTC\nM","MDD\nF","MDD\nM","BPD\nF","BPD\nM"))


p1 <- ggplot(demo, aes(x=cond_sex, y=RIN, color=condition))+
  ggbeeswarm::geom_beeswarm(cex=3, method = "center")+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent", width=.7)+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  labs(title="RIN", x="", y="RIN")+
  theme_minimal()+theme(legend.position="none", panel.grid.minor=element_blank())

p2 <- ggplot(demo, aes(x=cond_sex, y=BMI, color=condition))+
  ggbeeswarm::geom_beeswarm(cex=3, method = "center")+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent", width=.7)+
  ylim(0,80)+
  labs(y="kg/m2", x="", title="BMI")+
  theme_minimal()+theme(legend.position="none", panel.grid.minor=element_blank())


demo2 = group_by(demo, cond_sex) %>% add_tally(name="n_subjects") %>%
  group_by(cond_sex, condition, n_subjects, Smoking) %>% tally(name="n_smokers") %>%
  filter(Smoking=="Yes") %>% mutate(prop_smokers = n_smokers/n_subjects)

p3 <- ggplot(demo2, aes(x=cond_sex, y=prop_smokers, fill=condition))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=cpList$dx.pal)+
  ylim(0,1)+
  labs(y="Prop. of subjects w/ history", x="", title="Smoking")+
  theme_minimal()+theme(legend.position="none", panel.grid.minor=element_blank())

p4 <- ggplot(demo2, aes(x=cond_sex, y=n_smokers, fill=condition))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=cpList$dx.pal)+
  ylim(0,20)+
  labs(y="# of subjects w/ history", x="", title="Smoking")+
  theme_minimal()+theme(legend.position="none", panel.grid.minor=element_blank())


pdf(file="plots/publication/supp_covariate-selection/demographics.pdf", height=3, width=3)
p1
p2
p3
p4
dev.off()

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
