setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})
cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

comparisons2 = comparisons[c(1,3,5,2,4,6)]

results_set="seurat-pc30"
la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         coef=factor(coef, levels=comparisons2),
         is_deg=gene_id %in% la.degs$gene_id & adj.P.Val<.05) 
summary(lat$logFC)
summary(-log10(lat$adj.P.Val))

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$transfer.bright)[c(1:4,8,5:7)],
		labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2),
         is_deg=gene_id %in% lr.degs$gene_id & adj.P.Val<.05) 
summary(lrt$logFC)
summary(-log10(lrt$adj.P.Val))

mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv") %>%
	group_by(gene_name) %>% slice_max(MeanRatio, n=1) %>%
	mutate(cellType.target= ifelse(MeanRatio<1.5, "multi", cellType.target))

color.palette = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"=cpList$low.res.light[["L3"]],
            cpList$low.res.light[c("Astro")], cpList$low.res.bright["Oligo"],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="darkgrey",
	"NS"="lightgrey"
)

plot.df = left_join(lat, mratio.sn[,c("gene_name","MeanRatio","cellType.target")], by="gene_name")
plot.df[is.na(plot.df$cellType.target),"cellType.target"] = "NS"
plot.df[!plot.df$is_deg, "cellType.target"] = "NS"

p1 <- ggplot(mutate(plot.df, cluster="L-A"), aes(x=logFC, y=-log10(adj.P.Val), color=cellType.target))+
  rasterize(geom_point(size=.1, shape=4), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, is_deg==T) %>% mutate(cluster="L-A"), size=.5), dpi=300)+
#  scale_color_manual(values=c("grey","black"), guide="none")+
  scale_color_manual(values=color.palette)+
  facet_grid(cols=vars(coef), rows=vars(cluster))+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        #aspect.ratio=1, 
	text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
	axis.title.x=element_blank(), plot.margin=margin(0,5.5,0,5.5,"pt"))


# solo version for main figures
plot.df2 = mutate(plot.df, sex.group=factor(coef, levels=comparisons))
p1.vert <- ggplot(plot.df2, aes(x=logFC, y=-log10(adj.P.Val), color=cellType.target))+
  geom_point(size=.1, shape=4)+
  geom_point(data=filter(plot.df2, is_deg==T), size=.5)+
#  scale_color_manual(values=c("grey","black"), guide="none")+
  scale_color_manual(values=color.palette)+
  facet_wrap(vars(sex.group), ncol=2)+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, legend.position="none",
        text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
        axis.title.x=element_blank())
ggsave(file="plots/publication/volcano_plots/main_seurat.pdf", p1.vert,
	height=3, width=2)

plot.df2 = left_join(lrt, mratio.sn[,c("gene_name","MeanRatio","cellType.target")], by="gene_name")
plot.df2[is.na(plot.df2$cellType.target),"cellType.target"] = "NS"
plot.df2[!plot.df2$is_deg, "cellType.target"] = "NS"

p2 <- ggplot(plot.df2, aes(x=logFC, y=-log10(adj.P.Val), color=cellType.target))+
  rasterize(geom_point(size=.1, shape=4), dpi=300)+
  rasterize(geom_point(data=filter(plot.df2, is_deg==T), size=.5), dpi=300)+
#  scale_color_manual(values=c("grey","black"), guide="none")+
  scale_color_manual(values=color.palette)+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,7))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3))

lay_mat= rbind(c(1,1,1,1,1,1), matrix(2, ncol=6, nrow=8))
ggsave(file=paste0("plots/publication/volcano_plots/", results_set, "_all-volcanoes_mean-ratio-colored.pdf"), 
       grid.arrange(p1+theme(legend.position="none"), p2+theme(legend.position="none"), layout_matrix=lay_mat),
       width=6, height=9)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()


