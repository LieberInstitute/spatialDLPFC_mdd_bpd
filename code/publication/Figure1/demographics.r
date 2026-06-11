setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

demo = read.csv("processed-data/publication/demographics.csv")

demo$cond_sex = factor(paste(demo$condition, demo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BD F","BD M"),
	labels=gsub(" ", "\n", c("NTC F","NTC M","MDD F","MDD M","BD F","BD M")))
demo$condition = factor(demo$condition, levels=c("NTC","MDD","BD"), labels=c("NTC","MDD","BPD"))

p1 <- ggplot(demo, aes(x=cond_sex, y=age, color=condition))+
  ggbeeswarm::geom_beeswarm(cex=3)+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent", width=.7)+
  ylim(0,75)+labs(title="Age (donor)", x="", y="years")+
  theme_minimal()+theme(legend.position="none", panel.grid.minor=element_blank())

p2 <- ggplot(demo, aes(x=cond_sex, y=PMI, color=condition, shape=sex))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent")+
  ylim(0,60)+labs(title="PMI", x="", y="hours")+
  theme_minimal()+theme(legend.position="none")


ggsave(file="plots/publication/Figure1/demographics_age.pdf", 
       #gridExtra::grid.arrange(p1, p2, ncol=1),
	p1,
       bg="white", height=3, width=3)

stop("Age boxplot only")
median(demo$age)
summary(demo$age)
mean(demo$age)
sd(demo$age)
group_by(demo, sex) %>% summarise(med_age=median(age), avg_age=mean(age))
group_by(demo, cond_sex) %>% summarise(med_age=median(age))
