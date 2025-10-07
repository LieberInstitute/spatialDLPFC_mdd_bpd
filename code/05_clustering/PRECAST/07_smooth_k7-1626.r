setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(ggspavis)
})

set.seed(123)
source("code/05_clustering/PRECAST/07-supp_smoother.r")

cpList = readRDS("plots/colorPalettes.rds")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")
spe$precast_k9_1626_f = ifelse(is.na(spe$precast_k9_1626), "missing", as.character(spe$precast_k9_1626))
spe$precast_k7_1626_f = ifelse(is.na(spe$precast_k7_1626), "missing", as.character(spe$precast_k7_1626))

#annotate precast clusters
source("code/05_clustering/PRECAST/04-supp_cluster-names.r")
for(i in c("precast_k9_1626", "precast_k7_1626")) {
  colData(spe)[[paste0(i,"_f")]] = factor(colData(spe)[[paste0(i,"_f")]], levels= c(annotations[[i]],"missing"), 
                                        labels= c(names(annotations[[i]]), "missing"))
}
table(spe$precast_k7_1626_f, useNA="ifany")

cdata = as.data.frame(colData(spe))
# smooth precast_k9_1663 for all samples
all.samples = unique(cdata$sample_id)
cdataList <- lapply(all.samples, function(x){
  load(paste0("processed-data/04_feature_selection/per-sample_spe-conservative/",x,".Rdata"))
  cdata_sub = cdata[colnames(tmp),]
  tmp$precast_k7_1626 = cdata_sub$precast_k7_1626_f
  
  #first pass is very aggressive for GABA and low UMI
  new = smoother(labels_curr=cdata_sub$precast_k7_1626_f,
                 locs=spatialCoords(tmp),
                 props=c("L1"=.8, "L2"=.8, "L3/4"=.8, "L5"=.8,"L6"=.8, "WM"=.5, "low UMI"=.05, "missing"=.05),
                 k=10,
                 max_iter=2)
  #second pass takes care of the stragglers
  new2 = smoother(labels_curr=new,
                  locs=spatialCoords(tmp),
                  props=c("L1"=.4, "L2"=.5, "L3/4"=.5, "L5"=.5,"L6"=.5, "WM"=.4, "low UMI"=.05),
                  k=6,
                  max_iter=2)
  tmp$precast_k7_1626 = new2
  return(as.data.frame(colData(tmp)))
})

# extract colData
new_cdata = do.call(rbind, cdataList)
cat("\n\n")
table(new_cdata$precast_k7_1626, useNA="ifany")

stopifnot(identical(rownames(cdata), rownames(new_cdata)))

cdata$smoothed_k7_1626 = new_cdata$precast_k7_1626

write.csv(cdata[,c("sample_id","brnum","sex","condition","slide","key",
                   "precast_k9_1626","precast_k7_1626","precast_k7_1626_f","smoothed_k7_1626")], 
          "processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names = T)
cat("\ncolData csv with cluster results saved to: processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv\n")


# plot smoothed results

stopifnot(identical(rownames(cdata), colnames(spe)))
spe$smoothed_k7_1626_f = factor(cdata$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))

#create spe subsets for plotting
#ideal samples
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

mod_spatialCoords1 = spatialCoords(spe_sub)
for (i in unique(spe_sub$sample_id)) {
  tmp = mod_spatialCoords1[colData(spe_sub)$sample_id==i,]
  mod_spatialCoords1[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords1[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

#poorer-quality slides
lowq.samples = c("V13B23-282","V13B23-301","V13B23-302","V13B23-310","V13B23-311","V13B23-403")
spe_sub2 = spe[,spe$slide %in% lowq.samples]
cat("\nDim low quality sample spe:\n")
dim(spe_sub2)

spe_sub2$facet_row = factor(spe_sub2$slide, levels= lowq.samples, labels= gsub("-","\n",lowq.samples))
spe_sub2$facet_col = factor(spe_sub2$array, levels= c("A1","B1","C1","D1"))

mod_spatialCoords2 = spatialCoords(spe_sub2)
for (i in unique(spe_sub2$sample_id)) {
  tmp = mod_spatialCoords2[colData(spe_sub2)$sample_id==i,]
  mod_spatialCoords2[colData(spe_sub2)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords2[colData(spe_sub2)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

save.name = "k7-1626-smooth"
fill.palette = c(cpList$smoothed.bright, "low UMI"="lightgrey")
p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords1[,1], y_coord=mod_spatialCoords1[,2],
                sample_id="sample_id", annotate="smoothed_k7_1626_f", point_size=.1)+
  scale_color_manual(paste0("PRECAST\n", save.name), values=fill.palette)+
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


p2 <- plotSpots(spe_sub2, x_coord=mod_spatialCoords2[,1], y_coord=mod_spatialCoords2[,2],
                sample_id="sample_id", annotate="smoothed_k7_1626_f", point_size=.1)+
  scale_color_manual(paste0("PRECAST\n", save.name), values=fill.palette)+
  facet_grid(rows=vars(facet_row), cols=vars(facet_col), switch="y")+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

pdf(file=paste0("plots/05_clustering/PRECAST/example-spot-plots_conservative_", save.name, ".pdf"), height=10, width=8)
p1
p2
dev.off()
cat("\nSpot plots saved to:",paste0("plots/05_clustering/PRECAST/example-spot-plots_conservative_", save.name, ".pdf"),"\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
