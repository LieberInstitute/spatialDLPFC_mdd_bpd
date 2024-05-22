args=commandArgs(TRUE)

library(dplyr)
library(ggplot2)
library(SpatialExperiment)
library(scater)
library(here)
load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

title.key= c("sum_3MAD.outlier_slide"="UMI counts (3MAD), block per slide",
"sum_3MAD.outlier_sample"="UMI counts (3MAD), block per sample",
"genes_3MAD.outlier_slide"="n genes (3MAD), block per slide",
"genes_3MAD.outlier_sample"="n genes (3MAD), block per sample",
"chrM.ratio_3MAD.outlier_slide"="mito % (3MAD), block per slide",
"chrM.ratio_3MAD.outlier_sample"="mito % (3MAD), block per sample",
"sum_local.outlier"="UMI counts (local outliers)",
"genes_local.outlier"="n genes (local outliers)",
"chrM.ratio_local.outlier"="mito % (local outliers)")

title.key = title.key[args[[1]]]

plot.var.key= c("sum_3MAD.outlier_slide"="lg10.sum",
"sum_3MAD.outlier_sample"="lg10.sum",
"genes_3MAD.outlier_slide"="lg10.genes",
"genes_3MAD.outlier_sample"="lg10.genes",
"chrM.ratio_3MAD.outlier_slide"="expr_chrM_ratio",
"chrM.ratio_3MAD.outlier_sample"="expr_chrM_ratio",
"sum_local.outlier"="lg10.sum",
"genes_local.outlier"="lg10.genes",
"chrM.ratio_local.outlier"="expr_chrM_ratio")

plot.var.key = plot.var.key[args[[1]]]

pdf.key = c("sum_3MAD.outlier_slide"="umi-counts-slide_3MAD",
"sum_3MAD.outlier_sample"="umi-counts-sample_3MAD",
"genes_3MAD.outlier_slide"="n-genes-slide_3MAD",
"genes_3MAD.outlier_sample"="n-genes-sample_3MAD",
"chrM.ratio_3MAD.outlier_slide"="mito-perc-slide_3MAD",
"chrM.ratio_3MAD.outlier_sample"="mito-perc-sample-3MAD",
"sum_local.outlier"="umi-counts_local",
"genes_local.outlier"="n-genes_local",
"chrM.ratio_local.outlier"="mito-perc_local")

pdf.key = pdf.key[args[[1]]]

d1 = as.data.frame(colData(spe)[,c("key", "slide","position", "brain", args[[1]])])
colnames(d1)[ncol(d1)] = "var1"
d2 = group_by(d1, slide, position, brain) %>% summarise(var1 = sum(var1))

p1 <- ggplot(d2, aes(x=slide, y=var1, color=position))+
  geom_text(aes(label=brain))+ #change label to brain
  scale_color_brewer(palette="Dark2")+
  labs(y="n outliers",title=title.key)+
  theme_bw()

l1 = unique(spe$slide)
names(l1) = l1
slide.list = lapply(l1, function(x) colData(spe)$slide==x)
violin.list = lapply(seq_along(slide.list), function(x) {
  sub.spe = spe[,slide.list[[x]]]
  plotColData(sub.spe, x="brain",
              y=plot.var.key, color_by=args[[1]], point_size=.5)+
    scale_color_manual(values=c("darkgrey","red"))+
    coord_cartesian(ylim=c(0,round(max(colData(spe)[,plot.var.key])*1.05,1)))+
    labs(x="",color="outlier", title=paste("Slide",names(slide.list)[x]))+
    theme_bw()
})
p2 = PRECAST::drawFigs(violin.list, layout.dim = c(3, 2), common.legend = TRUE, legend.position = "right", align = "hv")

pdf(here("plots", "03_QC", paste0("outliers_distribution-plot_",pdf.key,".pdf")), width=8, height=12)
gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(1,3))
dev.off()
