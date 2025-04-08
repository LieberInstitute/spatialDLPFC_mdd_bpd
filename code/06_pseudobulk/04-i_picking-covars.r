#set a sense for these data results
enrich = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-age-sex-slide.csv", row.names=1)
enrich2 = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-sex.csv", row.names=1)
enrich3 = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-age-detected-PMI-RIN-sex.csv", row.names=1)
enrich4 = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-age-ncells-PMI-RIN-sex.csv", row.names=1)
enrich5 = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-detected-sex-slide.csv", row.names=1)
enrich6 = read.csv("processed-data/06_pseudobulk/results_layer-enrichment_covars-ncells-sex-slide.csv", row.names=1)

head(enrich)
head(enrich2)
identical(enrich$gene, enrich2$gene)
par(mfrow=c(1,1))
plot(enrich$t_stat_Vasc, enrich2$t_stat_Vasc, xlab="covars = age, sex, slide", ylab="covars = sex", main="t_stat_Vasc")
plot(enrich$t_stat_Vasc, enrich5$t_stat_Vasc, xlab="covars = age, sex, slide", ylab="covars = detected, sex, slide", main="t_stat_Vasc")
plot(enrich$t_stat_Vasc, enrich3$t_stat_Vasc, xlab="covars = age, sex, slide", ylab="covars = age, detected, PMI, RIN, sex", main="t_stat_Vasc")
plot(enrich2$t_stat_Vasc, enrich3$t_stat_Vasc, xlab="covars = sex", ylab="covars = age, detected, PMI, RIN, sex", main="t_stat_Vasc")
plot(enrich4$t_stat_Vasc, enrich3$t_stat_Vasc, xlab="covars = age, ncells, PMI, RIN, sex", ylab="covars = age, detected, PMI, RIN, sex", main="t_stat_Vasc")

plot(enrich$fdr_Vasc, enrich2$fdr_Vasc, xlab="covars = age, sex, slide", ylab="covars = sex", main="fdr_Vasc")
plot(enrich$fdr_Vasc, enrich5$fdr_Vasc, xlab="covars = age, sex, slide", ylab="covars = detected, sex, slide", main="fdr_Vasc")
plot(enrich2$fdr_Vasc, enrich3$fdr_Vasc, xlab="covars = sex", ylab="covars = age, detected, PMI, RIN, sex", main="fdr_Vasc")

plot(enrich$logFC_Vasc, enrich2$logFC_Vasc, xlab="covars = age, sex, slide", ylab="covars = sex", main="logFC_Vasc")


library(dplyr)
filter(enrich, fdr_Vasc<.05) %>% nrow() #7143
filter(enrich2, fdr_Vasc<.05) %>% nrow() #7809

filter(layer_modeling_results$enrichment, fdr_Layer1<.0001, logFC_Layer1>1) %>% nrow() #458
filter(enrich, fdr_L1<.0001, t_stat_L1>5) %>% nrow() #367
filter(enrich2, fdr_L1<.0001, t_stat_L1>5) %>% nrow() #366
filter(enrich3, fdr_L1<.0001, t_stat_L1>5) %>% nrow() #420


filter(enrich, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #458
filter(enrich2, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #441

#### **** this is where you start and orient yourself [look at correlation of fdr values]
############ WHEN CONSIDERING DETECTED AND SEX, THERE IS LITTLE EFFECT OF ADDING LOW IMPACT COVARS LIKE AGE, PMI, AND RIN
############ WHEN CONSIDERING DETECTED AND SEX, THERE IS LITTLE EFFECT OF ADDING SLIDE, LIKELY BECAUSE SAMPLE ID CORR ACCOUNTS FOR MOST OF THIS EFFECT
filter(enrich3, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #491 #covars = age, detected, PMI, RIN, sex
filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #493 #covars = detected, sex, slide
length(intersect(filter(enrich3, fdr_Vasc<.0001, logFC_Vasc>1)$gene, filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1)$gene)) #490
length(union(filter(enrich3, fdr_Vasc<.0001, logFC_Vasc>1)$gene, filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1)$gene)) #494
490/494 #.992
plot(enrich3$fdr_Vasc, enrich5$fdr_Vasc, xlab="covars = age, detected, pmi, rin, sex", ylab="covars = detected, sex, slide", main="fdr_Vasc")

plot(enrich6$fdr_Vasc, enrich5$fdr_Vasc, xlab="covars = ncells, sex, slide", ylab="covars = detected, sex, slide", main="fdr_Vasc")
filter(enrich6, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #533 #covars = ncells, sex, slide
filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1) %>% nrow() #493 #covars = detected, sex, slide
length(intersect(filter(enrich6, fdr_Vasc<.0001, logFC_Vasc>1)$gene, filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1)$gene)) #486
length(union(filter(enrich6, fdr_Vasc<.0001, logFC_Vasc>1)$gene, filter(enrich5, fdr_Vasc<.0001, logFC_Vasc>1)$gene)) #540
486/540 #.9
