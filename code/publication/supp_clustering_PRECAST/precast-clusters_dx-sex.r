setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(gridExtra)
  library(ggrastr)
})


set.seed(123)
setAutoBlockSize(1e9)

cpList = readRDS("plots/colorPalettes.rds")
fill.palette = cpList$smoothed.bright

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM"),
	labels=c("L1","L2","L3.4","L5","L6","WM"))


table(spe$smoothed_k9_1663, useNA="ifany")


# proportion of annotations in samples (bar)
cdata = group_by(as.data.frame(colData(spe)), sample_id) %>% add_tally(name="n_total") %>%
  group_by(sample_id, condition, sex, n_total, smoothed_k9_1663) %>%
  tally(name="nspots") %>%
  mutate(prop_spots=nspots/n_total, 
         condition=factor(condition, levels=c("NTC","MDD","BPD")),
         cond_sex= factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
  )


tmp1 = group_by(cdata, smoothed_k9_1663) %>% mutate(r1=rank(prop_spots)) %>%
  select(smoothed_k9_1663, sample_id, r1)
tmp2 = tidyr::pivot_wider(tmp1, names_from="smoothed_k9_1663", values_from="r1", values_fill=0)
order1 = arrange(tmp2, WM, L6)
cdata = mutate(cdata, x_lab= factor(sample_id, levels=order1$sample_id))

blist <- lapply(c("NTC","MDD","BPD"), function(x) {
  tmp = filter(cdata, condition==x)
  ggplot(tmp,
         aes(x=x_lab, y=prop_spots, fill=smoothed_k9_1663))+
    rasterize(geom_bar(stat="identity", position="fill", color="black", linewidth=.1), dpi=300)+
    scale_y_continuous(expand=expansion(0), breaks=c(0,.5,1), labels=c("0",".5","1"))+
    scale_x_discrete("capture area", expand=expansion(0))+
    facet_wrap(vars(cond_sex), ncol=2, scales="free_x")+
    scale_fill_manual(values=fill.palette, guide="none")+
    theme_bw()+theme(axis.text.x= element_blank(), axis.title.y=element_blank(),
                     strip.background = element_rect(fill="white", color=NA),
                     panel.grid=element_blank(),
                     axis.ticks = element_line(color="black", linewidth=.3))
})


# six samples
vistoseg.samples = c("332-A1","308-D1",
  "023-D1","309-D1",
  "382-B1","382-A1")
vistoseg.samples = paste0("V13B23-", gsub("-","_", vistoseg.samples))
vistoseg.samples[[3]] = "V13Y10-023_D1"

spe_sub = spe[,spe$sample_id %in% vistoseg.samples]
spe_sub$cond_sex = factor(paste(spe_sub$condition, spe_sub$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
table(spe_sub$cond_sex)
spe_list <- lapply(c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"), function(x) {
	spe_sub[,spe_sub$cond_sex==x]
})

plist <- lapply(spe_list, function(x) {
  p = make_escheR(x) %>%
    add_fill(var="smoothed_k9_1663", point_size = .8)
  p+scale_fill_manual(values=cpList$smoothed.bright, guide="none")+
    labs(title=unique(x$cond_sex))+theme(plot.title=element_text(size=14, hjust=.5))
})


lmat = rbind(c(1,2),c(3,3))
pdf(file="plots/publication/Figure1/precast-clusters_dx-sex.pdf", width=3, height=4)
grid.arrange(rasterize(plist[[1]], dpi=300), rasterize(plist[[2]], dpi=300), 
             blist[[1]]+theme(legend.position="none"), layout_matrix=lmat)
grid.arrange(rasterize(plist[[3]], dpi=300), rasterize(plist[[4]], dpi=300), 
             blist[[2]]+theme(legend.position="none"), layout_matrix=lmat)
grid.arrange(rasterize(plist[[5]], dpi=300), rasterize(plist[[6]], dpi=300), 
             blist[[3]]+theme(legend.position="none"), layout_matrix=lmat)
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
