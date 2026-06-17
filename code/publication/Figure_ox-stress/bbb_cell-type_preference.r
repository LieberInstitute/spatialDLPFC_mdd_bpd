library(dplyr)
library(ggplot2)

mod_subset = c("GFAP","GLUL","MT1M","HSPA1A","IFITM3","A2M","CD74")

sn.col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
#            "Micro"=cpList$low.res.light[["L3"]],
            "Micro"="#C28658",
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
        "multi"="grey"
)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
tmp = filter(refined.modules, TF %in% mod_subset) %>%
	bind_rows(data.frame(TF=mod_subset, target=mod_subset, importance=1))

mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv") %>% 
	group_by(gene, gene_name) %>% slice_max(MeanRatio, n=1) %>% 
	mutate(cellType.target= ifelse(MeanRatio<1.5, "multi", cellType.target),
		cellType.target = factor(cellType.target, levels=rev(c("Vasc","Astro","Oligo","Micro","InhN","ExcN","multi"))))


tmp2 <- left_join(tmp, mratio.sn[,c("gene","gene_name","cellType.target")], by=c("target"="gene_name")) %>% 
	filter(!is.na(cellType.target))

tmp3 = group_by(tmp2, TF, cellType.target) %>% tally() %>%
	mutate(TF=factor(TF, levels=mod_subset))

p1 <- ggplot(tmp3, aes(x=TF, y=n, fill=cellType.target))+
	geom_bar(stat="identity", position="stack", color="black", linewidth=.3)+
	scale_fill_manual(values=sn.col.pal)+
	theme_minimal()+theme(text=element_text(size=6))
	#coord_polar(theta="y")+
	#facet_wrap(vars(TF), scales="free")+
	#theme_void()+theme(aspect.ratio=1)

ggsave(file="plots/publication/Figure4/bbb_mean-ratio_bar-plot.pdf", p1, width=3, height=2)
