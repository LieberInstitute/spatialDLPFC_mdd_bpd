setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(dplyr)
        library(ggplot2)
	library(ggspavis)
})

set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
seurat_pc = "pc20"

cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$precast_k9_1663 = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI"),
	labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI"))

#low.res.pal = c("Micro/Vasc"="#911223","Astro"="#cfa45c",
#	"L2"="#5D9940", "L3"="#5095CD", 
#	"L4"="#85A0A0",
#        "L5"="#ddc94e","L6"="#E45C5F",
#        "Oligo"="#D1C4B0",
#        "Inhb"="#9377AC")

res = read.csv(paste0("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-",
	seurat_pc, "_red-precast-kweight-50-low-res.csv"), row.names=1)

# remove 10 extra spots with precast smoothed labels
#cdata = cdata[rownames(res),]
stopifnot(identical(rownames(cdata), rownames(res)))

cdata$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	########## merging pc30 L2 and L3
	#labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))
	labels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))

cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), 
                        levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))

######### merging pc30 L2 and L3
#seurat_pc = "pc30"
#low.res.pal = c("Micro/Vasc"="#911223","Astro"="#cfa45c",
#        "L2/3"="#5D9940", 
#	#"L3"="#5095CD",
#        "L4"="#85A0A0",
#        "L5"="#ddc94e","L6"="#E45C5F",
#        "Oligo"="#D1C4B0",
#        "Inhb"="#9377AC")

#bar plot of spot abundance
p1 <- ggplot(cdata, aes(x=cond_sex, fill=seurat_label))+
  geom_bar()+scale_fill_manual(paste0("Seurat\nlabel\n(",seurat_pc,")"), values=cpList$low.res.bright)+
  labs(y="# of spots", title=paste("Cell type abundance:",seurat_pc))+
  theme_minimal()+theme(axis.title.x=element_blank(), text=element_text(size=10), aspect.ratio=1, 
                        legend.position="bottom", plot.margin = margin(1.5,.5,1,.5, "cm"),
                        legend.key.size= unit(10,"pt"))


#plot heatmap
cdata2 = group_by(cdata, precast_k9_1663, seurat_label, .drop=F) %>% tally() %>%
  group_by(precast_k9_1663) %>% mutate(total_n=sum(n)) %>%
  ungroup() %>%
  mutate(prop_n=n/total_n,
         precast = factor(precast_k9_1663, levels=rev(levels(cdata$precast_k9_1663))))

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
plist <- lapply(names(cpList$low.res.bright), function(x) {
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


#plot low UMI proportions
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex),
                        levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
p1 <- ggplot(filter(cdata, precast_k9_1663=="low UMI"), aes(x=sample_id, fill=seurat_label))+
  geom_bar(stat="count", position="stack", color="black", linewidth=.3)+
  facet_wrap(vars(cond_sex), scales="free_x", ncol=2)+
  scale_fill_manual(values=cpList$low.res.bright)+
  labs(title="Only PRECAST k=9 low UMI cluster spots", fill=seurat_pc)+
  theme_bw()+theme(axis.text.x=element_blank())

plist = lapply(levels(cdata$seurat_label), function(x) {
	ggplot(filter(cdata, seurat_label==x), aes(x=sample_id, fill=precast_k9_1663))+
		geom_bar(stat="count", position="stack", color="black", linewidth=.3)+
		facet_wrap(vars(cond_sex), scales="free_x", ncol=2)+
		scale_fill_manual("PRECAST",values=c(cpList$earthy.pal2, "low UMI"="grey50"))+
		labs(title=paste("Only Seurat", seurat_pc, x))+
		theme_bw()+theme(axis.text.x=element_blank())
})

plist = c(list(p1), plist)
ggsave(file=paste0("plots/05_clustering/Seurat/per-sample-nspots_low-UMI-cluster_", seurat_pc, ".pdf"),
	gridExtra::marrangeGrob(grobs=plist, nrow=1, ncol=1, top=NULL),
	width=7.5, height=9)
cat("\n\nLow UMI assignemtn bar plot saved to:",paste0("plots/05_clustering/Seurat/per-sample-nspots_low-UMI-cluster_", seurat_pc,".png"),"\n")

print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

