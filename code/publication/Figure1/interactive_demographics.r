library(dplyr)
library(ggplot2)

cpList <- readRDS("plots/colorPalettes.rds")
source("code/02_build_spe/getMBvSampleInfo_function.r")

fastqs = grep("^V", list.files('raw-data/FASTQ'), value=T)
space = grep("^V", list.files("processed-data/01_spaceranger"), value=T)

demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")
head(demo)
#remove sample missing spacerange (run with NAc i think)
demo = demo[demo$sample_id!="V13B23-339_A1",]

demo$condition = factor(demo$condition, levels=c("NTC","MDD","BPD"))
demo$sex = factor(demo$sex, levels=c("F","M"))

#current sample list at 128 because of re-run slides
demo = demo[!demo$slide %in% c("V13Y10-020","V13B23-331"),]

dim(demo) #120

#remove repeated NTC M 
demo = demo[!demo$sample_id=="V13F27-338_C1",]

demo$cond_sex = factor(paste(demo$condition, demo$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))

p1 <- ggplot(demo, aes(x=cond_sex, y=age, color=condition, shape=sex))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent", width=.8)+
  ylim(0,75)+labs(title="Age (donor)", x="", y="years")+
  theme_minimal()+theme(legend.position="none")

p2 <- ggplot(demo, aes(x=cond_sex, y=PMI, color=condition, shape=sex))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$dx.pal)+
  geom_boxplot(outliers=F, color="black", fill="transparent")+
  ylim(0,60)+labs(title="PMI", x="", y="hours")+
  theme_minimal()+theme(legend.position="none")


ggsave(file="plots/publication/demographics.png", 
       gridExtra::grid.arrange(p1, p2, ncol=1),
       bg="white", height=8, width=4)

median(demo$age)
summary(demo$age)
mean(demo$age)
sd(demo$age)
group_by(demo, sex) %>% summarise(med_age=median(age), avg_age=mean(age))
group_by(demo, cond_sex) %>% summarise(med_age=median(age))
