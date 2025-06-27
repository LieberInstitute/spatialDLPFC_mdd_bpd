library(SpatialExperiment)
library(edgeR)

set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
spe_sub = spe_pseudo[,spe_pseudo$combined_cluster %in% c("L2","L3","L5")]
spe_sub$slide_name = gsub("-","\\.", spe_sub$slide)
spe_sub$combined_cluster = droplevels(spe_sub$combined_cluster)

d0 <- DGEList(counts(spe_sub))
d0 <- calcNormFactors(d0)

mm1 <- model.matrix(~0 + condition + combined_cluster + sex, colData(spe_sub))
is.fullrank(mm1)

#test without voom (skipping duplicateCorrelation for time/ testing purposes)
fit0 <- lmFit(logcounts(spe_sub), mm1)
head(coef(fit0))

cont_mtx = matrix(0, ncol=3, nrow=ncol(coef(fit0)), dimnames=list(colnames(coef(fit0)), c("NTC.MDD","NTC.BPD","MDD.BPD")))
cont_mtx["conditionNTC",c("NTC.MDD","NTC.BPD")] <- c(-1, -1)
cont_mtx["conditionMDD",c("NTC.MDD","MDD.BPD")] <- c(1, -1)
cont_mtx["conditionBPD",c("NTC.BPD","MDD.BPD")] <- c(1, 1)
cont_mtx

tmp0 = contrasts.fit(fit0, cont_mtx)
is.fullrank(tmp0$design)

eb0 <- eBayes(tmp0)
hist(topTable(eb0, "NTC.MDD", n=Inf)$P.Value) #good
hist(topTable(eb0, "NTC.BPD", n=Inf)$P.Value) #good
hist(topTable(eb0, "MDD.BPD", n=Inf)$P.Value) #good
table(topTable(eb0, "NTC.MDD", n=Inf)$adj.P.Val<.05) #311
table(topTable(eb0, "NTC.MDD", n=Inf)$adj.P.Val<.01) #79

mdd.sig.genes0 = rownames(topTable(eb0, "NTC.MDD", p.value = .01, n = Inf))

#test with voom
y1 <- voom(d0, mm1, plot=T)
#corfit = duplicateCorrelation(y1, block=spe_sub$sample_id)

fit1 = lmFit(y1)
head(coef(fit1))

tmp1 = contrasts.fit(fit1, cont_mtx)
is.fullrank(tmp1$design)

eb1 <- eBayes(tmp1)
hist(topTable(eb1, "NTC.MDD", n=Inf)$P.Value) #good
hist(topTable(eb1, "NTC.BPD", n=Inf)$P.Value) #good
hist(topTable(eb1, "MDD.BPD", n=Inf)$P.Value) #good
table(topTable(eb1, "NTC.MDD", n=Inf)$adj.P.Val<.05) #431
table(topTable(eb1, "NTC.MDD", n=Inf)$adj.P.Val<.01) #93


mdd.sig.genes1 = rownames(topTable(eb1, "NTC.MDD", p.value = .01, n = Inf))

plot(topTable(eb0, "NTC.MDD", n=Inf)[union(mdd.sig.genes0, mdd.sig.genes1),"t"],
     topTable(eb1, "NTC.MDD", n=Inf)[union(mdd.sig.genes0, mdd.sig.genes1),"t"],
     xlab="covars= cluster, sex", ylab="voom; covars= cluster, sex")

#add ncells covariate
mm2 <- model.matrix(~0 + condition + combined_cluster + sex + ncells, colData(spe_sub))
is.fullrank(mm2)
y2 <- voom(d0, mm2, plot=T)
#corfit = duplicateCorrelation(y1, block=spe_sub$sample_id)

fit2 = lmFit(y2)
head(coef(fit2))

cont_mtx2 = matrix(0, ncol=3, nrow=ncol(coef(fit2)), dimnames=list(colnames(coef(fit2)), c("NTC.MDD","NTC.BPD","MDD.BPD")))
cont_mtx2["conditionNTC",c("NTC.MDD","NTC.BPD")] <- c(-1, -1)
cont_mtx2["conditionMDD",c("NTC.MDD","MDD.BPD")] <- c(1, -1)
cont_mtx2["conditionBPD",c("NTC.BPD","MDD.BPD")] <- c(1, 1)
cont_mtx2
tmp2 = contrasts.fit(fit2, cont_mtx2)
is.fullrank(tmp2$design)

eb2 <- eBayes(tmp2)
hist(topTable(eb2, "NTC.MDD", n=Inf)$P.Value) #good
hist(topTable(eb2, "NTC.BPD", n=Inf)$P.Value) #good
hist(topTable(eb2, "MDD.BPD", n=Inf)$P.Value) #good
table(topTable(eb2, "NTC.MDD", n=Inf)$adj.P.Val<.05) #384
table(topTable(eb2, "NTC.MDD", n=Inf)$adj.P.Val<.01) #384

mdd.sig.genes2 = rownames(topTable(eb2, "NTC.MDD", p.value = .01, n = Inf))

plot(topTable(eb1, "NTC.MDD", n=Inf)[union(mdd.sig.genes1, mdd.sig.genes2),"t"],
     topTable(eb2, "NTC.MDD", n=Inf)[union(mdd.sig.genes1, mdd.sig.genes2),"t"],
     xlab="voom; covars= cluster, sex", ylab="voom; covars= cluster, sex, ncells")


#add detected covariate
mm3 <- model.matrix(~0 + condition + combined_cluster + sex + ncells + detected, colData(spe_sub))
is.fullrank(mm3) #FALSE
y3 <- voom(d0, mm3, plot=T)

fit3 = lmFit(y3)
head(coef(fit3))

cont_mtx3 = matrix(0, ncol=3, nrow=ncol(coef(fit3)), dimnames=list(colnames(coef(fit3)), c("NTC.MDD","NTC.BPD","MDD.BPD")))
cont_mtx3["conditionNTC",c("NTC.MDD","NTC.BPD")] <- c(-1, -1)
cont_mtx3["conditionMDD",c("NTC.MDD","MDD.BPD")] <- c(1, -1)
cont_mtx3["conditionBPD",c("NTC.BPD","MDD.BPD")] <- c(1, 1)
cont_mtx3
tmp3 = contrasts.fit(fit3, cont_mtx3)
is.fullrank(tmp3$design) #FALSE

eb3 <- eBayes(tmp3)
hist(topTable(eb3, "NTC.MDD", n=Inf)$P.Value) #bad
hist(topTable(eb3, "NTC.BPD", n=Inf)$P.Value) #good
hist(topTable(eb3, "MDD.BPD", n=Inf)$P.Value) #bad? yes bad
table(topTable(eb3, "NTC.MDD", n=Inf)$adj.P.Val<.05) #1


#add detected covariate with no voom
fit4 = lmFit(logcounts(spe_sub), mm3)
head(coef(fit4))

tmp4 = contrasts.fit(fit4, cont_mtx3)
is.fullrank(tmp4$design) #FALSE

eb4 <- eBayes(tmp4)
hist(topTable(eb4, "NTC.MDD", n=Inf)$P.Value) #good
hist(topTable(eb4, "NTC.BPD", n=Inf)$P.Value) #good
hist(topTable(eb4, "MDD.BPD", n=Inf)$P.Value) #good
table(topTable(eb4, "NTC.MDD", n=Inf)$adj.P.Val<.05) #376
table(topTable(eb4, "NTC.MDD", n=Inf)$adj.P.Val<.01) #101

mdd.sig.genes4 = rownames(topTable(eb4, "NTC.MDD", p.value = .01, n = Inf))

plot(topTable(eb0, "NTC.MDD", n=Inf)[union(mdd.sig.genes0, mdd.sig.genes4),"t"],
     topTable(eb4, "NTC.MDD", n=Inf)[union(mdd.sig.genes0, mdd.sig.genes4),"t"],
     xlab="covars= cluster, sex", ylab="covars= cluster, sex, ncells, detected")
