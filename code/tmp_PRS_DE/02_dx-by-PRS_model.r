setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
})
set.seed(123)


prs_wrapper <- function(risk_set, include_RIN=T) {
# cutoffs
cutoffs = c("MDD"="1e.07", "Bipolar"="1e.06")

# load PRS
prs.df = read.csv(paste0("raw-data/PRS/PRS_",risk_set,".csv"))
colnames(prs.df)[1] = "brnum"
prs.df$risk = risk_set

prs.df = prs.df[,c("brnum","risk",paste0("p.cutoff.",cutoffs[[risk_set]]))]
colnames(prs.df)[3] = "PRS"
prs.df$PRS_scaled = scale(prs.df$PRS)

# load in donor metadata
cdata = read.csv("raw-data/sample_info/DLPFC_cross-disorders_demographics_MBv.csv") 
cdata$RIN_scaled = scale(cdata$RIN)

prs.df = left_join(prs.df, cdata) %>%
  mutate(sex=factor(sex, levels=c("F","M")), 
         condition=factor(condition, levels=c("NTC","MDD","BPD")))

# pull SNP PCs
t1 = read.table("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/tqtl_in/astro.gene.covars.txt", row.names=1)

t2 = as.data.frame(apply(t(t1[-1,]), MARGIN=2, as.numeric))
t2$brnum = as.character(t1[1,])

snp.df = t2[,c("brnum",grep("snp", colnames(t2), value=T))]
snp.df2 = apply(snp.df[,-1], MARGIN=2, scale)
colnames(snp.df2) = paste0(colnames(snp.df2),"_scaled")
snp.df = cbind(snp.df, snp.df2)
stopifnot(nrow(snp.df)==119)

# merge snp covars
all.df = left_join(prs.df, snp.df)

list1 = lapply(c("MDD","BPD"), function(i) {
	set.seed(123)
	cond.df = filter(all.df, condition %in% c("NTC",i)) %>% mutate(condition= droplevels(condition))
	if(include_RIN) m1 = glm(condition ~ PRS_scaled + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled + RIN_scaled, data=cond.df, family=binomial)
	if(!include_RIN) m1 = glm(condition ~ PRS_scaled + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled, data=cond.df, family=binomial)
	return(m1)
})
names(list1) = paste0("dx", c("MDD","BPD"))

all.df2 = mutate(all.df, condition= factor(condition, levels=c("NTC","MDD","BPD"), labels=c("NTC","DX","DX")))
if(include_RIN) list1[["dxAny"]] = glm(condition ~ PRS_scaled + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled + RIN_scaled, data=all.df2, family=binomial)
if(!include_RIN) list1[["dxAny"]] = glm(condition ~ PRS_scaled + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled, data=all.df2, family=binomial)

return(list1)
}

# include RIN in covars
cat("\n\nRUNNING VERSION WITH RIN INCLUDED AS COVARIATE...\n")
riskList = list("prsMDD"="MDD", "prsBPD"="Bipolar")
riskList = lapply(riskList, prs_wrapper)
saveRDS(riskList, "processed-data/tmp_PRS_DE/glm-binomial_MDD-Bipolar-PRS_predicting-dxMDD-dxBPD_with-RIN.rda")
cat("\nList of model output saved to: processed-data/tmp_PRS_DE/glm-binomial_MDD-Bipolar-PRS_predicting-dxMDD-dxBPD_with-RIN.rda\n")

cat("\n\ndxMDD ~ prsMDD + snpPCs + RIN\n")
print(summary(riskList[["prsMDD"]]$dxMDD))

cat("\n\ndxBPD ~ prsMDD + snpPCs + RIN\n")
print(summary(riskList[["prsMDD"]]$dxBPD))

cat("\n\ndxAny ~ prsMDD + snpPCs + RIN\n")
print(summary(riskList[["prsMDD"]]$dxAny))

cat("\n\ndxMDD ~ prsBPD + snpPCs + RIN\n")
print(summary(riskList[["prsBPD"]]$dxMDD))

cat("\n\ndxBPD ~ prsBPD + snpPCs + RIN\n")
print(summary(riskList[["prsBPD"]]$dxBPD))

cat("\n\ndxAny ~ prsBPD + snpPCs + RIN\n")
print(summary(riskList[["prsBPD"]]$dxAny))

# don't include RIN in covars
cat("\n\nRUNNING VERSION WITHOUT RIN INCLUDED AS COVARIATE...\n")
riskList = list("prsMDD"="MDD", "prsBPD"="Bipolar")
riskList = lapply(riskList, prs_wrapper, include_RIN=F)
saveRDS(riskList, "processed-data/tmp_PRS_DE/glm-binomial_MDD-Bipolar-PRS_predicting-dxMDD-dxBPD_without-RIN.rda")
cat("\nList of model output saved to: processed-data/tmp_PRS_DE/glm-binomial_MDD-Bipolar-PRS_predicting-dxMDD-dxBPD_without-RIN.rda\n")

cat("\n\ndxMDD ~ prsMDD + snpPCs\n")
print(summary(riskList[["prsMDD"]]$dxMDD))

cat("\n\ndxBPD ~ prsMDD + snpPCs\n")
print(summary(riskList[["prsMDD"]]$dxBPD))

cat("\n\ndxAny ~ prsMDD + snpPCs\n")
print(summary(riskList[["prsMDD"]]$dxAny))

cat("\n\ndxMDD ~ prsBPD + snpPCs\n")
print(summary(riskList[["prsBPD"]]$dxMDD))

cat("\n\ndxBPD ~ prsBPD + snpPCs\n")
print(summary(riskList[["prsBPD"]]$dxBPD))

cat("\n\ndxAny ~ prsBPD + snpPCs\n")
print(summary(riskList[["prsBPD"]]$dxAny))

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
