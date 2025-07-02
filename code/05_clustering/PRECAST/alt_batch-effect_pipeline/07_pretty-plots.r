setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(ggspavis)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
#spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
#stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))

#spe$precast_k9_1663 = factor(spotdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM","low UMI"))
#spe$precast_k7_1663_f = factor(spe$precast_k7_1663, 
#	levels=c(4,3,2,6,1,5,7),
#	labels=c("Vasc","Vasc/L1","L2","L3/4","L5","L6","WM"))

spe$precast_k7_HM_f = factor(spe$precast_k7_HM,
        levels=c(3,1,2,5,7,4,6),
        labels=c("Vasc","L1","L2","L3/4","L5","L6","WM"))

vistoseg.samples = c("308-D1", "329-A1", "342-B1", "332-B1", #NTC M
                     "342-A1", "332-A1", "327-C1", "329-B1",#NTC F
                     "382-C1", "309-D1","352-A1", "329-C1",  #MDD M
                     "382-D1", "352-B1", "309-C1", "380-B1", #"329-D1", #MDD F
                     "382-A1", "352-C1", "342-D1", "327-B1", #BPD M
                     "382-B1", "308-A1", "309-A1", "327-A1") #BPD F

vistoseg.samples = paste0("V13B23-", gsub("-","_", vistoseg.samples))

spe_sub = spe[,spe$sample_id %in% vistoseg.samples]
cat("\nDim VistoSeg sample spe:\n")
dim(spe_sub)

spe_sub$cond_sex = factor(paste(spe_sub$condition, spe_sub$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))
spe_sub$facet_col = factor(spe_sub$sample_id, levels=vistoseg.samples, labels=paste0("C",rep(c(1:4), 6)))

mod_spatialCoords = spatialCoords(spe_sub)
for (i in unique(spe_sub$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

cpList = readRDS("plots/colorPalettes.rds")
#fill.palette = c(cpList[["earthy.pal2"]],"low UMI"="grey90")
#fill.palette = c("Vasc"=cpList$earthy.pal2[["Vasc"]], "Vasc/L1"=cpList$earthy.pal2[["L1"]],
#	"L2"=cpList$earthy.pal2[["L2"]], "L3/4"=cpList$earthy.pal2[["L3"]],
#	cpList$earthy.pal2[c("L5","L6","WM")])
fill.palette = c("Vasc"=cpList$earthy.pal2[["Vasc"]], "L1"=cpList$earthy.pal2[["L1"]],
        "L2"=cpList$earthy.pal2[["L2"]], "L3/4"=cpList$earthy.pal2[["L3"]],
        cpList$earthy.pal2[c("L5","L6","WM")])
p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
                sample_id="sample_id", annotate="precast_k7_HM_f", point_size=.1)+
  scale_color_manual("PRECAST\ncluster", values=fill.palette)+
  facet_grid(rows=vars(cond_sex), cols=vars(facet_col), switch="y")+
  geom_text(aes(x=0, y= -60, 
                label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13,
                              .data[["sample_id"]], "")), hjust=0, vjust=0, 
            color="black", size=2)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

#ggsave("plots/05_clustering/PRECAST/VistoSeg-examples_n1663-k9_spot-plots.png", 
#       p1,
#       bg="white", width=8, height=10)
#cat("\nSpot plots of high quality, representative sections from all dx and sex groups saved to: plots/05_clustering/PRECAST/VistoSeg-examples_n1663-k9_spot-plots.png\n")

#poorer-quality slides
lowq.samples = c("V13B23-282","V13B23-301","V13B23-302","V13B23-310","V13B23-311","V13B23-403")
spe_sub2 = spe[,spe$slide %in% lowq.samples]
cat("\nDim low quality sample spe:\n")
dim(spe_sub2)

spe_sub2$facet_row = factor(spe_sub2$slide, levels= lowq.samples, labels= gsub("-","\n",lowq.samples))
spe_sub2$facet_col = factor(spe_sub2$array, levels= c("A1","B1","C1","D1"))

mod_spatialCoords = spatialCoords(spe_sub2)
for (i in unique(spe_sub2$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub2)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub2)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub2)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

p2 <- plotSpots(spe_sub2, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
                sample_id="sample_id", annotate="precast_k7_HM_f", point_size=.1)+
  scale_color_manual("PRECAST\ncluster", values=fill.palette)+
  facet_grid(rows=vars(facet_row), cols=vars(facet_col), switch="y")+
  #geom_text(aes(x=0, y= -60, 
  #              label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13,
  #                            .data[["sample_id"]], "")), hjust=0, vjust=0, 
  #          color="black", size=2)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

#ggsave("plots/05_clustering/PRECAST/low-qual-examples_n1663-k9_spot-plots.png", 
#       p2,
#       bg="white", width=8, height=10)
#cat("\nSpot plots of lower quality slides with high # low UMI cluster saved to: plots/05_clustering/PRECAST/low-qual-examples_n1663-k9_spot-plots.png\n")

pdf(file="plots/05_clustering/PRECAST/example-spot-plots_HM-k7.pdf", height=10, width=8)
p1
p2
dev.off()
cat("\nSpot plots saved to: plots/05_clustering/PRECAST/example-spot-plots_n1663-k7.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
