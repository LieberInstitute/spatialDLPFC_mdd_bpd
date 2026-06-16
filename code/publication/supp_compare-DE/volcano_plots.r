setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons


source("code/09_DEG_GRN/load_DEGs.r")
comparisons2 = comparisons[c(1,3,5,2,4,6)]


# domain-SP
plot.df = filter(de.df, source=="sm") %>% mutate(sex.group= factor(sex.group, levels=comparisons2))

p1 <- ggplot(filter(plot.df, cluster=="L-A"), aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(size=.1, color="grey"), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, cluster=="L-A", adj.P.Val2<.05), size=.1, color="black"), dpi=300)+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        #aspect.ratio=1, 
	text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
	axis.title.x=element_blank(), plot.margin=margin(0,5.5,0,5.5,"pt"))


p2 <- ggplot(filter(plot.df, cluster!="L-A"), aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(size=.1, color="grey"), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, cluster!="L-A", adj.P.Val2<.05), size=.1, color="black"), dpi=300)+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,7))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3))

lay_mat= rbind(c(1,1,1,1,1,1), matrix(2, ncol=6, nrow=6))
ggsave(file="plots/publication/supp_compare-DE/domain-SP_all-volcanoes.pdf", 
       grid.arrange(p1+theme(legend.position="none"), p2+theme(legend.position="none"), layout_matrix=lay_mat),
       width=6, height=7)


# domain-CT
plot.df = filter(de.df, source=="se") %>% mutate(sex.group= factor(sex.group, levels=comparisons2),
	cluster=factor(cluster, levels=c("L-A","Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
		labels=c("L-A","M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))) 

p1 <- ggplot(filter(plot.df, cluster=="L-A"), aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(size=.1, color="grey"), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, cluster=="L-A", adj.P.Val2<.05), size=.1, color="black"), dpi=300)+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        #aspect.ratio=1, 
        text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
        axis.title.x=element_blank(), plot.margin=margin(0,5.5,0,5.5,"pt"))

p2 <- ggplot(filter(plot.df, cluster!="L-A"), aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(size=.1, color="grey"), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, cluster!="L-A", adj.P.Val2<.05), size=.1, color="black"), dpi=300)+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,7))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3))


lay_mat= rbind(c(1,1,1,1,1,1), matrix(2, ncol=6, nrow=8))
ggsave(file="plots/publication/supp_compare-DE/domain-CT_all-volcanoes.pdf",
       grid.arrange(p1+theme(legend.position="none"), p2+theme(legend.position="none"), layout_matrix=lay_mat),
       width=6, height=9)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
