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
seurat_pc = "pc30"

cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$precast_k9_1663 = factor(cdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI"),
	labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI"))


res = read.csv(paste0("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-",
	seurat_pc, "_red-precast-kweight-50-low-res.csv"), row.names=1)

# remove 10 extra spots with precast smoothed labels
#cdata = cdata[rownames(res),]
stopifnot(identical(rownames(cdata), rownames(res)))

cdata$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	########## merging pc30 L2 and L3
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))

cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), 
                        levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))

p1 <- ggplot(filter(cdata, precast_k9_1663=="low UMI"), aes(x=sample_id, fill=seurat_label))+
  geom_bar(stat="count", position="stack", color="black", linewidth=.1)+
  facet_wrap(vars(cond_sex), scales="free_x", ncol=2)+
  scale_y_continuous(expand=expansion(0))+
  scale_fill_manual(values=cpList$transfer.bright, guide="none")+
  labs(fill=seurat_pc, y="# spots")+
  theme_minimal()+theme(axis.text.x=element_blank(), text=element_text(size=8),
	panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
	axis.ticks = element_line(color="black", linewidth=.3))

ggsave(file="plots/publication/supp_clustering_Seurat/seurat-label_low-UMI-proportion.pdf", p1, height=4, width=3)


cdata2 <- group_by(cdata, sample_id, condition, sex) %>% add_tally(name="n_total") %>%
  group_by(sample_id, condition, sex, n_total, seurat_label) %>%
  tally(name="nspots") %>%
  mutate(prop_spots=nspots/n_total, 
         condition=factor(condition, levels=c("NTC","MDD","BPD")),
         cond_sex= factor(paste(condition, sex), levels=c("NTC F","MDD F","BPD F","NTC M","MDD M","BPD M"))
  )


tmp1 = group_by(cdata2, seurat_label) %>% mutate(r1=rank(prop_spots)) %>%
  select(seurat_label, sample_id, r1)
tmp2 = tidyr::pivot_wider(tmp1, names_from="seurat_label", values_from="r1", values_fill=0)
order1 = arrange(tmp2, Oligo, L6)
cdata2 = mutate(cdata2, x_lab= factor(sample_id, levels=order1$sample_id))

p2 <- ggplot(cdata2, aes(x=x_lab, y=prop_spots, fill=seurat_label))+
    geom_bar(stat="identity", position="fill", color="black", linewidth=.1)+
    scale_y_continuous(expand=expansion(0), breaks=c(0,.5,1), labels=c("0",".5","1"))+
    scale_x_discrete("capture area", expand=expansion(0))+
    facet_wrap(vars(cond_sex), ncol=3, scales="free_x")+
    scale_fill_manual(values=cpList$transfer.bright, guide="none")+
    theme_bw()+theme(axis.text.x= element_blank(), axis.title.y=element_blank(),
                     strip.background = element_rect(fill="white", color=NA),
                     panel.grid=element_blank(), text=element_text(size=8),
                     axis.ticks = element_line(color="black", linewidth=.3))

ggsave(file="plots/publication/supp_clustering_Seurat/seurat-label_per-sample-proportion.pdf", p2, height=2, width=3.5)

print("\n\nReproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
