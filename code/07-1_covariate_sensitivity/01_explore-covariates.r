setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

results_set = "smoothed-n1663-k9"
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("smoothed", colnames(colData(spe_pseudo)))] = "cluster"

#results_set = "seurat-pc30"
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
#colnames(colData(spe_pseudo))[grep("seurat", colnames(colData(spe_pseudo)))] = "cluster"

cat("\nResults set:", results_set, "\n")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]

cdata = read.csv("raw-data/sample_info/Page_BMI_Smoking_CrossDisorder_120825.csv")
#colSums(is.na(cdata))

cat("\n\nNew donor data filtered to MBv donors\n\n")
cdata = filter(cdata, BrNum %in% spe_pseudo$brnum)
table(cdata[,c("Nicotine","Cotinine")])
#.............Cotinine
#Nicotine     Negative Not Tested Positive
#..Negative         54          0        8 #8 negative for nicotine and positive for cotinine
#..Not Tested        0         37        0
#..Positive          0          0       20
cat("\n\n")
table(cdata[,c("Smoking","Cotinine")])
#.............Cotinine
#Nicotine     Negative Not Tested Positive
#..Negative         54          0        8
#..Not Tested        0         37        0
#..Positive          0          0       20

#decide to use Smoking status
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
#dim(sdata)

sdata = left_join(sdata, cdata[,c("BrNum","BMI","Smoking")], by=c("brnum"="BrNum"))
#colSums(is.na(sdata))

#co-linearity? not an issue
cat("\n\nCorrelation between continuous donor covariates:\n")
cor(sdata[,c("age","PMI","RIN","BMI")])
#.....age         PMI         RIN         BMI
#age  1.00000000  0.17342499  0.05166935 -0.07249217
#PMI  0.17342499  1.00000000  0.06303647 -0.16666976
#RIN  0.05166935  0.06303647  1.00000000 -0.16601974
#BMI -0.07249217 -0.16666976 -0.16601974  1.00000000

cat("\n\nSmoking status is extremely confounded with diagnosis:\n")
table(sdata[,c("condition","sex","Smoking")])
#, , Smoking = No
#...........sex
#condition  F  M
#......NTC 17 17
#......MDD  7  7
#......BPD  4  8

#, , Smoking = Yes
#...........sex
#condition  F  M
#......NTC  3  3
#......MDD 13 12
#......BPD 16 12

#plot if needed
if(!file.exists("plots/07-1_covariate_sensitivity/donor-covariates_by_dx-sex.png")) {
o1 <- ggplot(sdata, aes(x=sex, y=age, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width = .75)+
  geom_boxplot(fill="white", alpha=.5)+
  scale_color_manual(values=cpList$dx.pal)+
  ylim(0,65)+labs(y="years", title="Age")+theme_bw()

o2 <- ggplot(sdata, aes(x=sex, y=BMI, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width = .75)+
  geom_boxplot(fill="white", alpha=.5)+
  scale_color_manual(values=cpList$dx.pal)+
  ylim(0,80)+
  labs(y="kg/m2", title="BMI")+theme_bw()

o3 <- ggplot(sdata, aes(x=sex, y=RIN, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width = .75)+
  geom_boxplot(fill="white", alpha=.5)+
  scale_color_manual(values=cpList$dx.pal)+
  geom_hline(aes(yintercept=7), lty=2)+
  ylim(5,10)+
  labs(y="RIN", title="RIN")+theme_bw()

o4 <- ggplot(sdata, aes(x=sex, y=PMI, color=condition))+
  ggbeeswarm::geom_beeswarm(dodge.width = .75)+
  geom_boxplot(fill="white", alpha=.5)+
  scale_color_manual(values=cpList$dx.pal)+
  ylim(0,60)+
  labs(y="hours", title="PMI")+theme_bw()

sdata2 = group_by(sdata, sex, condition) %>% add_tally(name="n_subjects") %>%
  group_by(sex, condition, n_subjects, Smoking) %>% tally(name="n_smokers") %>%
  filter(Smoking=="Yes") %>% mutate(prop_smokers = n_smokers/n_subjects)

o5 = ggplot(sdata2, aes(x=sex, y=prop_smokers, fill=condition))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=cpList$dx.pal)+
  ylim(0,1)+
  labs(y="Prop. of subjects w/ history", title="Smoking")+theme_bw()

ggsave(file="plots/07-1_covariate_sensitivity/donor-covariates_by_dx-sex.png",
	grid.arrange(o1, o2, o3, o4, o5, ncol=2),
	width=7, height=9, bg="white")
cat("\nDonor covariate boxplots saved to: plots/07-1_covariate_sensitivity/donor-covariates_by_dx-sex.png\n")
}


# update colData
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
stopifnot(identical(spe_pseudo$total, new.cdata$total))
colData(spe_pseudo) <- new.cdata

# pca with new variables
bio.vars = c("cluster","sum","detected","nspots","chrM_ratio")
bio.colors = RColorBrewer::brewer.pal(n=length(bio.vars), "Set1")
names(bio.colors) = bio.vars

exp.vars = c("sample_id", "condition", "sex", "slide", "seq","pc3")
exp.colors = RColorBrewer::brewer.pal(n=length(exp.vars), "Dark2")
names(exp.colors) = exp.vars

donor.vars = c("age", "BMI", "Smoking", "RIN", "PMI")
donor.colors = c("#E4775D","#A0C255","#EEBC4A","grey","#93D3F6")
names(donor.colors) = donor.vars

p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=exp.vars)+
        scale_y_continuous()+
	scale_color_manual("", values= exp.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
	labs(title="Experimental variables", y="% PC variance explained")

p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=bio.vars)+
        scale_y_continuous()+
        scale_color_manual("", values=bio.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
        labs(title="Domain variables", y="% PC variance explained")

p3 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=10, variables=donor.vars)+
        scale_y_continuous()+
        scale_color_manual("", values=donor.colors)+
	scale_x_continuous(breaks=c(2,4,6,8,10))+
        labs(title="Donor variables", y="% PC variance explained")

#percent variance
var.m = getVarianceExplained(spe_pseudo, variables=c(exp.vars, bio.vars, donor.vars), exprs_values="logcounts")

# exp vars
tmp1 = tibble::rownames_to_column(as.data.frame(var.m[,exp.vars]), var="gene_id") %>%
	tidyr::pivot_longer(all_of(exp.vars), names_to="variable", values_to="variance")
tmp1.summ = group_by(tmp1, variable) %>% summarise(avg.var=mean(variance)) %>%
	arrange(desc(avg.var))
tmp1.summ$variable = factor(tmp1.summ$variable, levels=tmp1.summ$variable)
tmp1.summ$y_pos = seq(from=.75, by=-.15, length.out=nrow(tmp1.summ))
tmp1.summ$text = paste0(tmp1.summ$variable," = ", round(tmp1.summ$avg.var, 2), "%")

p4 <- ggplot(tmp1, aes(x=variance, color=variable))+
	stat_ecdf(linewidth=2)+
	geom_text(data=tmp1.summ, aes(x=95, y=y_pos, label=text, color=variable), hjust=1)+
	geom_text(data=data.frame("y_pos"=.90, "text"="Avg. gene var."), aes(x=95, y=y_pos, label=text), color="black", hjust=1)+
	scale_color_manual(values=exp.colors, guide="none")+
	xlim(0,100)+
	labs(x="% gene variance explained", title="Experimental variables", y="ecdf")+
	theme_minimal()+theme(plot.title=element_text(size=10))

#bio vars
tmp1 = tibble::rownames_to_column(as.data.frame(var.m[,bio.vars]), var="gene_id") %>%
        tidyr::pivot_longer(all_of(bio.vars), names_to="variable", values_to="variance")
tmp1.summ = group_by(tmp1, variable) %>% summarise(avg.var=mean(variance)) %>%
        arrange(desc(avg.var))
tmp1.summ$variable = factor(tmp1.summ$variable,	levels=tmp1.summ$variable) 
tmp1.summ$y_pos = seq(from=.75, by=-.15, length.out=nrow(tmp1.summ))
tmp1.summ$text = paste0(tmp1.summ$variable," = ", round(tmp1.summ$avg.var, 2), "%")

p5 <- ggplot(tmp1, aes(x=variance, color=variable))+
        stat_ecdf(linewidth=2)+
        geom_text(data=tmp1.summ, aes(x=95, y=y_pos, label=text, color=variable), hjust=1)+
        geom_text(data=data.frame("y_pos"=.90, "text"="Avg. gene var."), aes(x=95, y=y_pos, label=text), color="black", hjust=1)+
        scale_color_manual(values=bio.colors, guide="none")+
	xlim(0,100)+
        labs(x="% gene variance explained", title="Domain variables", y="ecdf")+
        theme_minimal()+theme(plot.title=element_text(size=10))

#donor vars
tmp1 = tibble::rownames_to_column(as.data.frame(var.m[,donor.vars]), var="gene_id") %>%
        tidyr::pivot_longer(all_of(donor.vars), names_to="variable", values_to="variance")
tmp1.summ = group_by(tmp1, variable) %>% summarise(avg.var=mean(variance)) %>%
        arrange(desc(avg.var))
tmp1.summ$variable = factor(tmp1.summ$variable,	levels=tmp1.summ$variable) 
tmp1.summ$y_pos = seq(from=.75, by=-.15, length.out=nrow(tmp1.summ))
tmp1.summ$text = paste0(tmp1.summ$variable," = ", round(tmp1.summ$avg.var, 2), "%")

p6 <- ggplot(tmp1, aes(x=variance, color=variable))+
        stat_ecdf(linewidth=2)+
        geom_text(data=tmp1.summ, aes(x=95, y=y_pos, label=text, color=variable), hjust=1)+
        geom_text(data=data.frame("y_pos"=.90, "text"="Avg. gene var."), aes(x=95, y=y_pos, label=text), color="black", hjust=1)+
        scale_color_manual(values=donor.colors, guide="none")+
	xlim(0,100)+
        labs(x="% gene variance explained", title="Donor variables", y="ecdf")+
        theme_minimal()+theme(plot.title=element_text(size=10))


#heatmap of gene variance
cor.var.m = cor(var.m, method="spearman")
hmp = pheatmap::pheatmap(cor.var.m, angle_col=90, silent=T, 
	main=paste0("Gene var. correlation (",results_set, ")"),
	treeheight_col=10)

lay_mat= rbind(c(1,4),c(2,5),c(3,6),c(7,7),c(7,7))
ggsave(file=paste0("plots/07-1_covariate_sensitivity/explore-covariates_", results_set,".png"),
	grid.arrange(p1, p2, p3, p4, p5, p6, hmp[[4]], layout_matrix=lay_mat),
	width=8, height=10, bg="white")
cat("\nExploratory plots saved to:", paste0("plots/07-1_covariate_sensitivity/explore-covariates_", results_set,".png"),"\n")


#subset by cluster and see variance explained correlation
varList <- lapply(levels(spe_pseudo$cluster), function(x) {
  scater::getVarianceExplained(spe_pseudo[,spe_pseudo$cluster==x], 
                               variables=c(exp.vars[-1],bio.vars[-1],donor.vars), 
                               exprs_values="logcounts")
})
names(varList) <- levels(spe_pseudo$cluster)

if(results_set=="smoothed-n1663-k9") clus.colors = cpList$smoothed.bright
if(results_set=="seurat-pc30") clus.colors = cpList$transfer.bright

corList <- lapply(varList, function(x) {
  cor(x, method="spearman")
})

hmpList = lapply(names(corList), function(x) {
  tmp = pheatmap::pheatmap(corList[[x]], angle_col = 90, main=x, silent=T,
                     treeheight_row = 15, treeheight_col = 15)
  return(tmp[[4]])
})

pdf(file=paste0("plots/07-1_covariate_sensitivity/explore-covariates_", results_set, "_by-cluster.pdf"), width=8, height=11)
par(mfrow=c(5,3))
for(i in c(exp.vars[-1],bio.vars[-1],donor.vars)) {
  for(j in names(varList)) {
    if(j==names(varList)[[1]]) {
      plot(ecdf(varList[[j]][,i]), col=clus.colors[[j]], main=i, xlim=c(0,100),
           xlab="% variance explained")
    } else {
      lines(ecdf(varList[[j]][,i]), col=clus.colors[[j]])
    }
  }
}
do.call(grid.arrange, c(grobs=hmpList, ncol=2))
dev.off()
cat("\nExploratory plots (by cluster) saved to:", paste0("plots/07-1_covariate_sensitivity/explore-covariates_", results_set,"_by-cluster.png"),"\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
