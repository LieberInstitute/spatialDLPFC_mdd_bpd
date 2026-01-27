library(dplyr)
library(ggplot2)

cpList <- readRDS("plots/colorPalettes.rds")

#prs scores: https://pmc.ncbi.nlm.nih.gov/articles/PMC7612115/

prs.df <- do.call(rbind, lapply(c("MDD","Bipolar","SCZ"), function(x) {
  tmp = read.csv(paste0("raw-data/PRS/PRS_",x,".csv"))
  colnames(tmp)[1] = "brnum"
  tmp$risk = x
  return(as.data.frame(tmp))
})) %>%
  mutate(risk=factor(risk, levels=c("MDD","Bipolar","SCZ"), labels=c("MDD","BPD","SCZ")))
head(prs.df)

prs.df2 = prs.df[,c(1:10,ncol(prs.df))]
prs.df2 = tidyr::pivot_longer(prs.df2, all_of(colnames(prs.df2)[2:10]),
                              names_to="p.cutoff", names_prefix = "p.cutoff.",
                              values_to="score") %>%
  mutate(p.cutoff= factor(p.cutoff, levels=gsub("p\\.cutoff\\.", "", colnames(prs.df2[2:10]))
                          )
         )

# load in donor metadata
cdata = read.csv("raw-data/sample_info/DLPFC_cross-disorders_demographics_MBv.csv")                 

prs.df2 = left_join(prs.df2, cdata) %>%
  mutate(sex=factor(sex, levels=c("F","M")), 
         condition=factor(condition, levels=c("NTC","MDD","BPD")))

p1 <- ggplot(filter(prs.df2, risk=="MDD"), aes(x=sex, y=score, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=cpList$dx.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: MDD")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))

p2 <- ggplot(filter(prs.df2, risk=="BPD"), aes(x=sex, y=score, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=cpList$dx.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: Bipolar")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))

p3 <- ggplot(filter(prs.df2, risk=="SCZ"), aes(x=sex, y=score, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=cpList$dx.pal)+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: SCZ")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))


pdf(file="plots/09_PRS_DE/PRS_p-cutoffs.pdf")
p1
p2
p3
dev.off()


bp.df = read.csv("raw-data/sample_info/Bipolar_subtype.csv")
head(bp.df)

prs.df2 = left_join(prs.df2, bp.df[,c("BrNum","BP.subtype")], by=c("brnum"="BrNum")) %>%
  mutate(condition2= factor(ifelse(is.na(BP.subtype), as.character(condition), BP.subtype), 
                            levels=c("NTC","MDD","BP1","BP2")))
head(prs.df2)

p1 <- ggplot(filter(prs.df2, risk=="MDD"), aes(x=condition2, y=score, color=condition2))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=c("NTC"="#7f7f7f", "MDD"="#FB8861FF", "BP1"="#9260b2", "BP2"="#291f43"))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: MDD")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))

p2 <- ggplot(filter(prs.df2, risk=="BPD"), aes(x=condition2, y=score, color=condition2))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=c("NTC"="#7f7f7f", "MDD"="#FB8861FF", "BP1"="#9260b2", "BP2"="#291f43"))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: Bipolar")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))

p3 <- ggplot(filter(prs.df2, risk=="SCZ"), aes(x=condition2, y=score, color=condition2))+
  ggbeeswarm::geom_beeswarm(dodge.width=.75)+
  geom_boxplot(fill="white", alpha=.5, outliers=F)+
  scale_color_manual(values=c("NTC"="#7f7f7f", "MDD"="#FB8861FF", "BP1"="#9260b2", "BP2"="#291f43"))+
  scale_y_continuous(labels=function(x) format(x, scientific=T, digits=2))+
  facet_wrap(vars(p.cutoff), ncol=3, scales="free_y")+
  labs(title="GWAS set: SCZ")+theme(aspect.ratio=1, axis.text.y=element_text(size=7))


pdf(file="plots/09_PRS_DE/PRS_p-cutoffs_bp-subtype.pdf")
p1
p2
p3
dev.off()



prs.df3 = bind_rows(filter(prs.df2, risk=="MDD", p.cutoff=="1e.07") %>% 
                      select(brnum, age, sex, condition, PMI, RIN, BP.subtype, risk, score),
                    filter(prs.df2, risk=="BPD", p.cutoff=="1e.06") %>% 
                      select(brnum, age, sex, condition, PMI, RIN, BP.subtype, risk, score),
                    filter(prs.df2, risk=="SCZ", p.cutoff=="1e.08") %>%
                      select(brnum, age, sex, condition, PMI, RIN, BP.subtype, risk, score)) %>%
  tidyr::pivot_wider(names_from="risk", values_from="score")

ggplot(prs.df3, aes(x=MDD, y=BPD, color=condition))+
  geom_point()+
  scale_color_manual(values=cpList$dx.pal)+
  facet_wrap(vars(sex))+
  scale_x_continuous("PRS (MDD)", labels=function(x) format(x, scientific=T, digits=2))+
  scale_y_continuous("PRS (BPD)", labels=function(x) format(x, scientific=T, digits=2))+
  theme(aspect.ratio=1)

colnames(prs.df3)[8:10] = c("MDD_p.1e07","BPD_p.1e06","SCZ_p.1e08")
write.csv(prs.df3, "raw-data/PRS/PRS_chosen-p-cutoffs.csv", row.names=F)
