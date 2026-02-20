library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrastr)

cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

comparisons2 = comparisons[c(1,3,5,2,4,6)]

results_set="smoothed-k9-1663"
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
p1 <- ggplot(lat, aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
  rasterize(geom_point(size=.1), dpi=300)+
  rasterize(geom_point(data=filter(lat, is_deg==T), size=.1), dpi=300)+
  scale_color_manual(values=c("grey","black"), guide="none")+
  facet_wrap(vars(coef), ncol=3)+
  coord_cartesian(xlim=c(-2.5,2.5))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6),
                        panel.border = element_rect(fill=NA, color="black", linewidth=.5))

ggsave(file="plots/publication/volcanoes.pdf", 
       p1,
       width=3, height=3)

lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                      results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
         cluster=factor(cluster, levels=names(cpList$smoothed.bright)),
         sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2),
         is_deg=gene_id %in% lr.degs$gene_id & adj.P.Val<.05) 
summary(lrt$logFC)
summary(-log10(lrt$adj.P.Val))
summary(filter(lrt, cluster=="L1")$logFC)
plist <- lapply(names(cpList$smoothed.bright), function(x) {
  ggplot(filter(lrt, cluster==x), aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
    rasterize(geom_point(size=.1), dpi=300)+
    rasterize(geom_point(data=filter(lrt, is_deg==T, cluster==x), size=.1), dpi=300)+
    scale_color_manual(values=c("grey","black"), guide="none")+
    facet_wrap(vars(sex.group), ncol=3)+
    coord_cartesian(xlim=c(-4,4), ylim=c(0,6))+
    theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                          aspect.ratio=1, text=element_text(size=6),
                          panel.border = element_rect(fill=NA, color=cpList$smoothed.bright[[x]], linewidth=.5))
})

plist[[7]] <- p1
names(plist) <- c(names(cpList$smoothed.bright),"L-A")
plist = plist[c(7,1:6)]

ggsave(file="plots/publication/Figure2/volcanoes.pdf", 
       marrangeGrob(grobs = plist, ncol=1, nrow=1, top = quote(names(plist)[g])),
       width=3, height=3)

p2 <- ggplot(lrt, aes(x=logFC, y=-log10(adj.P.Val), color=is_deg))+
  rasterize(geom_point(size=.1), dpi=300)+
  rasterize(geom_point(data=filter(lrt, is_deg==T), size=.1), dpi=300)+
  scale_color_manual(values=c("grey","black"), guide="none")+
  facet_grid(cols=vars(sex.group), rows=vars(cluster))+
  coord_cartesian(xlim=c(-4,4), ylim=c(0,6))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3))

ggsave(file="plots/publication/Figure2/volcanoes_L-R-supp.pdf", 
       p2,
       width=6, height=6)
