setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(dplyr)
        library(ggplot2)
	library(ggspavis)
})

set.seed(123)

seurat_pc = "pc30"

cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$precast_k9_1663_f = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","L5","L6","WM","GABA","low UMI"))

low.res.pal = c("Micro/Vasc"="#911223","Astro"="#cfa45c",
	"L2"="#5D9940", "L3"="#5095CD", 
	"L4"="#85A0A0",
        "L5"="#ddc94e","L6"="#E45C5F",
        "Oligo"="#D1C4B0",
        "Inhb"="#9377AC")

res = read.csv(paste0("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-",
	seurat_pc, "_red-precast-kweight-50-low-res.csv"), row.names=1)
stopifnot(identical(rownames(cdata), rownames(res)))
cdata$seurat_label = factor(res$predicted.id, levels=names(low.res.pal),
	########## merging pc30 L2 and L3
	labels=c("Micro/Vasc","Astro","L2/3","L2/3","L4","L5","L6","Oligo","Inhb"))
	##########
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), 
                        levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))

######### merging pc30 L2 and L3
seurat_pc = "pc30-L2-L3-merge"
low.res.pal = c("Micro/Vasc"="#911223","Astro"="#cfa45c",
        "L2/3"="#5D9940", 
	#"L3"="#5095CD",
        "L4"="#85A0A0",
        "L5"="#ddc94e","L6"="#E45C5F",
        "Oligo"="#D1C4B0",
        "Inhb"="#9377AC")

#bar plot of spot abundance
p1 <- ggplot(cdata, aes(x=cond_sex, fill=factor(seurat_label, levels=names(low.res.pal))))+
  geom_bar()+scale_fill_manual(paste0("Seurat\nlabel\n(",seurat_pc,")"),values=low.res.pal)+
  labs(y="# of spots", title=paste("Cell type abundance:",seurat_pc))+
  theme_minimal()+theme(axis.title.x=element_blank(), text=element_text(size=10), aspect.ratio=1, 
                        legend.position="bottom", plot.margin = margin(1.5,.5,1,.5, "cm"),
                        legend.key.size= unit(10,"pt"))


#plot heatmap
cdata2 = group_by(cdata, precast_k9_1663_f, seurat_label, .drop=F) %>% tally() %>%
  group_by(precast_k9_1663_f) %>% mutate(total_n=sum(n)) %>%
  ungroup() %>%
  mutate(prop_n=n/total_n,
         precast = factor(precast_k9_1663_f, levels=rev(levels(cdata$precast_k9_1663_f))))

p2 <- ggplot(cdata2, aes(x=seurat_label, y=precast, fill=log10(n+1)))+
  geom_tile(color="grey90")+
  scale_fill_gradientn("# spots\n(log10)",colors=colorRampPalette(c("white","grey70","black"), bias=.5)(6))+
  geom_text(data=mutate(cdata2, prop_n=round(prop_n,2)*100) %>% filter(prop_n>15),
            aes(label=prop_n), color="red3", size=3, fontface="bold")+
  scale_x_discrete(paste("Label transfer: seurat", seurat_pc), expand = c(0,0))+
  scale_y_discrete("Original PRECAST k=9 annotation", expand = c(0,0))+
  ggtitle(paste("MBv label transfer:", seurat_pc))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   text=element_text(size=10),
                   aspect.ratio=1, legend.key.size = unit(10, "pt"))

#plot individual distribution spot plots
load("processed-data/04_feature_selection/per-sample_spe/V13B23-329_A1.Rdata")
cdata_tmp = cdata[colnames(tmp),]
stopifnot(identical(rownames(cdata_tmp), rownames(colData(tmp))))
plist <- lapply(names(low.res.pal), function(x) {
  tmp$plot.me = cdata_tmp$seurat_label==x
  plotSpots(tmp, annotate="plot.me", point_size = .5)+
    scale_color_manual(values=c("grey","red3"), guide="none")+
    ggtitle(x)
})

#laymat = rbind(c(1,1,1,2,2,2),
#               c(1,1,1,2,2,2),
#               c(3,3,4,4,5,5),
#               c(6,6,7,7,8,8),
#               c(9,9,10,10,11,11))
########## merging pc30 L2 and L3
laymat = rbind(c(1,1,1,2,2,2),
               c(1,1,1,2,2,2),
               c(3,3,4,4,NA,NA),
               c(5,5,6,6,7,7),
               c(8,8,9,9,10,10))
########## 

ggsave(file=paste0("plots/05_clustering/Seurat/MBv_label-transfer-", seurat_pc,"_plots.png"), 
	gridExtra::grid.arrange(p1, p2, 
		plist[[1]], plist[[2]], plist[[3]], 
		plist[[4]], plist[[5]], plist[[6]], 
		plist[[7]], plist[[8]], 
		########## merging pc30 L2 and L3
		#plist[[9]], 
		##########
		layout_matrix=laymat),
	bg="white", width=10, height=14)
cat("\n\nResults plots saved to:",paste0("plots/05_clustering/Seurat/MBv_label-transfer-", seurat_pc,"_plots.png"),"\n")


#examples spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
stopifnot(identical(rownames(cdata), colnames(spe)))
spe$seurat_label = cdata$seurat_label
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
low.res.pal[["L4"]] = "#c2cfcf"
p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords1[,1], y_coord=mod_spatialCoords1[,2],
	sample_id="sample_id", annotate="seurat_label", point_size=.1)+
	scale_color_manual(paste0("Seurat\nlabel\n(", seurat_pc,")"), values=low.res.pal)+
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
	scale_color_manual(paste0("Seurat\nlabel\n(", seurat_pc,")"), values=low.res.pal)+
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
