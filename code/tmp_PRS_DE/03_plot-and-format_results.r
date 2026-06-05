args = commandArgs(TRUE)
cat("\nGWAS set for PRS predictor:", args[[1]],"\n")
x = paste0("prs", args[[1]])

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)
source("code/07_dx_DE/custom_functions.r")
cpList <- readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
results_set = "smoothed-k9-1663"

# add PRS
prs.group = c("prsMDD","prsBPD","prsSCZ")
cdata = read.csv("raw-data/PRS/PRS_chosen-p-cutoffs.csv")
cdata = cdata[,c(1:2, 5:6, 8:10)]
colnames(cdata)[5:7] = prs.group
cdata[,5:7] = scale(cdata[,5:7])

# join with other metadata needed for plotting
sdata = distinct(as.data.frame(colData(spe_pseudo)[,c("brnum","sample_id","condition","sex","age","PMI","RIN")]))
sdata = left_join(sdata, cdata, by=c("brnum","age","PMI","RIN"))

# PRS boxplot
sdata$condition = factor(sdata$condition, levels=c("NTC","MDD","BPD"))
sdata$plot_me = sdata[[x]]

p1 <- ggplot(sdata, aes(x=condition, y=plot_me))+
  ggbeeswarm::geom_beeswarm(aes(color=condition, shape=sex), cex=3)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+scale_shape_manual(values=c("F"=16, "M"=17))+
  geom_boxplot(outliers=F, fill="transparent")+ylim(-3,3)+
  labs(x="diagnosis", y="PRS (z-score)")+
  theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

p2 <- ggplot(mutate(sdata, x_group=factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"),
	labels=gsub(" ","\n",c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")))), 
    aes(x=x_group, y=plot_me, color=condition))+
  geom_boxplot(fill="white")+ylim(-3,3)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+
  labs(y="PRS (z-score)")+
  theme_bw()+theme(axis.title.x=element_blank(), 
	panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# merge PRS data with spe
new.cdata = merge(colData(spe_pseudo), sdata, sort=F)
stopifnot(identical(spe_pseudo$total, new.cdata$total)) 
colData(spe_pseudo) <- new.cdata

plist = list(p1, p2)
for (results_set in c("smoothed-k9-1663","seurat-pc30")) {
	# if needed, load in spe object with correct gene set gene names
	if(results_set=="seurat-pc30") load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
# load in PRS model results
results <- readRDS(paste0("processed-data/tmp_PRS_DE/lmFit-voom_layer-adjusted_", 
                          results_set, 
                          "_", x, "-sexM_rev-gene-input_covars-pc3-age-nspots.rda"))

# make contrast matrix
c.mx = cbind("F_x" = as.numeric(colnames(results$coefficients)==x),
             "M_x" = as.numeric(colnames(results$coefficients)==x))
rownames(c.mx) = colnames(results$coefficients)
c.mx[paste0(x,":sexM"),"M_x"] = 1

# eBayes fit
eb <- eBayes(contrasts.fit(results, c.mx), trend=T, robust=T)

# F test
eb.f = getTopTable(eb, .coef="all") %>% mutate("F_design"=paste0("~ ", x, "*sexM"))
f_sig = filter(eb.f, adj.P.Val<.05)$gene_id

# moderated t test
sex.res = sexTopTable(eb, phist=F) %>% mutate(group=x)
sex.res$adj.P.Val = p.adjust(sex.res$P.Value, method="BH")

cat("\nLayer-adjusted analysis with adj p<.05 and F adj p<.05:\n")
print(filter(sex.res, adj.P.Val<.05, gene_id %in% f_sig) %>% 
	group_by(sex, group, .drop=F) %>% tally())

phist = ggplot(sex.res, aes(x=P.Value))+
  geom_histogram(bins=50)+facet_grid(cols=vars(sex), rows=vars(group))+
  labs(title=paste0(results_set, ": T-test p val histogram"))+
  theme_bw()

p3 <- ggplot(sex.res, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3)+
  geom_point(data=filter(sex.res, adj.P.Val<.05, gene_id %in% f_sig), size=.3, color="red2")+
  facet_grid(rows=vars(group), cols=vars(sex))+
  xlim(-ceiling(max(sex.res$logFC)), ceiling(max(sex.res$logFC)))+
  ggtitle(paste0("Layer-adjusted (", results_set, ")"))+
  theme_bw()

plist = c(plist, list(phist, p3))

write.csv(eb.f, paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set,
	"_rev-gene-input_F-test.csv"), row.names=F)
cat("\n\n\nSaved F test results dframe to:", paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set,
        "_rev-gene-input_F-test.csv"), "\n")
write.csv(sex.res, paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set, 
	"_rev-gene-input_moderated-t-test.csv"), row.names=F)
cat("\nSaved moderated t test results dframe to:", paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set, 
	"_rev-gene-input_moderated-t-test.csv"), "\n\n")
}

#plist <- list(p1, p2, phist, p3)
lay_mat = rbind(c(1,2),c(3,3),c(4,4),c(5,5),c(6,6))
ggsave(file=paste0("plots/tmp_PRS_DE/", x, "_model-summaries.pdf"),
	arrangeGrob(grobs=plist, layout_matrix=lay_mat, top=x),
	width=6, height=10)

cat("\n\nSaved plots to:", paste0("plots/tmp_PRS_DE/", x, "_model-summaries.pdf"),"\n\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
