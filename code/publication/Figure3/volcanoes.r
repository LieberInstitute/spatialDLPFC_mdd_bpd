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
#p1 <- ggplot(lat, aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
#  rasterize(geom_point(size=.1), dpi=300)+
#  rasterize(geom_point(data=filter(lat, is_deg==T), size=.1), dpi=300)+
#  scale_color_manual(values=c("grey","black"), guide="none")+
#  facet_wrap(vars(coef), ncol=3)+
#  coord_cartesian(xlim=c(-2.5,2.5))+
#  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
#                        aspect.ratio=1, text=element_text(size=6),
#                        panel.border = element_rect(fill=NA, color="black", linewidth=.5))
#
#ggsave(file="plots/publication/volcanoes.pdf", 
#       p1,
#       width=3, height=3)

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$transfer.bright)[c(1:4,8,5:7)],
		labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2),
         is_deg=gene_id %in% lr.degs$gene_id & adj.P.Val<.05) 
summary(lrt$logFC)
summary(-log10(lrt$adj.P.Val))
#stop("Check L-A max adj. p value.")
#summary(filter(lrt, cluster=="L1")$logFC)
#plist <- lapply(names(cpList$smoothed.bright), function(x) {
#  ggplot(filter(lrt, cluster==x), aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
#    rasterize(geom_point(size=.1), dpi=300)+
#    rasterize(geom_point(data=filter(lrt, is_deg==T, cluster==x), size=.1), dpi=300)+
#    scale_color_manual(values=c("grey","black"), guide="none")+
#    facet_wrap(vars(sex.group), ncol=3)+
#    coord_cartesian(xlim=c(-4,4), ylim=c(0,6))+
#    theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
#                          aspect.ratio=1, text=element_text(size=6),
#                          panel.border = element_rect(fill=NA, color=cpList$smoothed.bright[[x]], linewidth=.5))
#})
#
#plist[[7]] <- p1
#names(plist) <- c(names(cpList$smoothed.bright),"L-A")
#plist = plist[c(7,1:6)]
#
#ggsave(file="plots/publication/Figure2/volcanoes.pdf", 
#       marrangeGrob(grobs = plist, ncol=1, nrow=1, top = quote(names(plist)[g])),
#       width=3, height=3)

p1 <- ggplot(mutate(lat, cluster="L-A"), aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
  rasterize(geom_point(size=.1), dpi=300)+
  rasterize(geom_point(data=filter(lat, is_deg==T) %>% mutate(cluster="L-A"), size=.1), dpi=300)+
  scale_color_manual(values=c("grey","black"), guide="none")+
  facet_grid(cols=vars(coef), rows=vars(cluster))+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        #aspect.ratio=1, 
	text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
	axis.title.x=element_blank(), plot.margin=margin(0,5.5,0,5.5,"pt"))


p2 <- ggplot(lrt, aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
  rasterize(geom_point(size=.1), dpi=300)+
  rasterize(geom_point(data=filter(lrt, is_deg==T), size=.1), dpi=300)+
  scale_color_manual(values=c("grey","black"), guide="none")+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,7))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3))

lay_mat= rbind(c(1,1,1,1,1,1), matrix(2, ncol=6, nrow=8))
ggsave(file="plots/publication/Figure3/supp_seurat-pc30_all-volcanoes.pdf", 
       grid.arrange(p1, p2, layout_matrix=lay_mat),
       width=6, height=9)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
