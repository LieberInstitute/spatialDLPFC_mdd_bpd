setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
        library(dplyr)
        library(ggplot2)
        library(ggspavis)
})
set.seed(123)
setAutoBlockSize(1e9)

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")


cpList <- readRDS("plots/colorPalettes.rds")
seurat_pc = "pc20"

res = read.csv(paste0("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-",
        seurat_pc, "_red-precast-kweight-50-low-res.csv"), row.names=1)
spe = spe[,rownames(res)]
res = res[colnames(spe),]
stopifnot(identical(rownames(colData(spe)), rownames(res)))

## remove low UMI cluster for L4 and Inhb
#cdata = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)
#cdata = cdata[rownames(res),]
#cdata$predicted.id = res$predicted.id
#remove.spots = cdata$smoothed_k7_1626=="low UMI" & cdata$predicted.id %in% c("L4","Inhb")
#cdata2 = cdata[!remove.spots,]
#
#spe = spe[,rownames(cdata2)]
#res = res[colnames(spe),]

spe$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
#        ########## merging pc30 L2 and L3
#        labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))
	labels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))
spe$cond_sex = factor(paste(spe$condition, spe$sex),
                        levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))

#seurat_pc = "pc30-L2-L3-merge"

# proportion of annotations in samples (bar)

cdata = group_by(as.data.frame(colData(spe)), sample_id) %>% add_tally(name="n_total") %>%
  group_by(sample_id, condition, sex, n_total, seurat_label) %>%
  tally(name="nspots") %>%
  mutate(prop_spots=nspots/n_total,
         condition=factor(condition, levels=c("NTC","MDD","BPD")),
         cond_sex= factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
  )


tmp1 = group_by(cdata, seurat_label) %>% mutate(r1=rank(prop_spots)) %>%
  select(seurat_label, sample_id, r1)
tmp2 = tidyr::pivot_wider(tmp1, names_from="seurat_label", values_from="r1", values_fill=0)
order1 = arrange(tmp2, Oligo, L6)

split.name = unlist(strsplit(seurat_pc, "-"))
p1 <- ggplot(mutate(cdata, x_lab= factor(sample_id, levels=order1$sample_id)),
       aes(x=x_lab, y=prop_spots, fill=seurat_label))+
  geom_bar(stat="identity", position="fill", color="black", linewidth=.3)+
  facet_wrap(vars(cond_sex), ncol=2, scales="free_x")+
  scale_fill_manual(values=cpList$low.res.bright)+
  labs(title="Seurat label", x="sample ID", y="prop. of spots",
        fill=paste0("Seurat\n", seurat_pc))+
	#fill=paste0("Seurat ", split.name[1], "\n", paste(split.name[2:3], collapse="/"), " ", split.name[4]))+
  theme_bw()+theme(axis.text.x= element_blank())

ggsave(file=paste0("plots/05_clustering/Seurat/per-sample-proportions_", seurat_pc, ".png"),
       p1,
       bg="white", width=7, height=9)



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


#generate and save plots
#split.name = unlist(strsplit(seurat_pc, "-"))
p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords1[,1], y_coord=mod_spatialCoords1[,2],
	sample_id="sample_id", annotate="seurat_label", point_size=.1)+
	scale_color_manual(paste0("Seurat ", seurat_pc),
		#paste0("Seurat ", split.name[1], "\n", paste(split.name[2:3], collapse="/"), " ", split.name[4]), 
		values=cpList$low.res.bright)+
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
	sample_id="sample_id", annotate="seurat_label", point_size=.1)+
	scale_color_manual(paste0("Seurat ", seurat_pc),
		#paste0("Seurat ", split.name[1], "\n", paste(split.name[2:3], collapse="/"), " ", split.name[4]), 
		values=cpList$low.res.bright)+
	facet_grid(rows=vars(facet_row), cols=vars(facet_col), switch="y")+
	theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
		strip.text.x = element_blank(),
		strip.text.y.left = element_text(angle=0),
		legend.title=element_text(margin=margin(0,0,4,0,"pt")),
		legend.box.spacing = unit(2,"pt"),
		legend.margin=margin(0,0,0,0,"pt"),
		legend.box.margin = margin(0,2,0,2,"pt"))

pdf(file=paste0("plots/05_clustering/Seurat/example-spot-plots_MBv_label-transfer-", seurat_pc, ".pdf"), height=10, width=8)
	print(p1)
	print(p2)
dev.off()
cat("\nSpot plots saved to:",paste0("plots/05_clustering/Seurat/example-spot-plots_MBv_label-transfer-", seurat_pc, ".pdf"),"\n")

print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
