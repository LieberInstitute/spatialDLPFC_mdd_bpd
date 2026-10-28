setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
        library(dplyr)
        library(ggplot2)
})

set.seed(123)


cpList <- readRDS("plots/colorPalettes.rds")

load("processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res_norm-filt.Rdata")
sce_mdd <- sce_pseudo

load("processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res_norm-filt.Rdata")
sce_bd <- sce_pseudo

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
spe_mbv <- spe_pseudo

# assemble full dataframe across 3 studies
full.df = bind_rows(select(as.data.frame(colData(sce_mdd)), sample_id=Sample, condition, sex, seurat_low.res,
			nobs=ncells, detected, sum, subsets_ribo_percent) %>% 
		mutate(source="Maitra", clus_orig= seurat_low.res,
		seurat_label = factor(seurat_low.res, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"),
		labels= c("M.V","Ast","L2.3","L2.3","L4","Inb","L5","L6","Olg"))),
	select(as.data.frame(colData(sce_bd)), sample_id=individualID, condition, sex, seurat_low.res,
			nobs=ncells, detected, sum, subsets_ribo_percent) %>% 
		mutate(source="SZBD", clus_orig= seurat_low.res,
		seurat_label = factor(seurat_low.res, levels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"),
		labels= c("M.V","Ast","L2.3","L2.3","L4","Inb","L5","L6","Olg"))),
	select(as.data.frame(colData(spe_mbv)), sample_id, condition, sex, seurat_label,
			nobs=nspots, detected, sum, subsets_ribo_percent) %>% 
		mutate(source="MBv",
		seurat_label = factor(seurat_label, levels=levels(spe_mbv$seurat_label),
			labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")),
		clus_orig= seurat_label,
		seurat_low.res = factor(seurat_label, levels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"),
			labels=c("M.V","Ast","L2","L4","Inb","L5","L6","Olg")))
) 

# to accommodate L2.3 merge in SRT, duplicate the L2.3 entries to appear in both L2 and L3 facets
full.df = bind_rows(full.df, 
		filter(full.df, source=="MBv", seurat_label=="L2.3") %>% mutate(seurat_low.res="L3", clus_orig=NA)) %>%
	mutate(clus_orig = factor(clus_orig, levels=c("M.V","Ast","L2","L2.3","L3","L4","Inb","L5","L6","Olg")),
		seurat_low.res = factor(seurat_low.res, levels=levels(sce_mdd$seurat_low.res)),
		lg10.sum= log10(sum),
		condition= factor(condition, levels=c("NTC","MDD","BD","BPD"), labels=c("NTC","MDD","BPD","BPD")),
		x_lab= factor(paste(source, condition, sep="\n"), levels=c("MBv\nNTC","MBv\nMDD","MBv\nBPD",
			"Maitra\nNTC","Maitra\nMDD","SZBD\nNTC","SZBD\nBPD")))
# easy labeling of this nuance in plots
text.df = data.frame(x_pos= c(1,1), y_pos= c(1750,1750), seurat_low.res=factor(c("L2","L3"), levels=levels(full.df$seurat_low.res)), tlabel=rep("MBv values\nfor L2/3", 2))

# plot QC metrics
p1 <- ggplot(full.df, aes(x=x_lab, y=nobs, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=cpList$dx.pal)+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res), scales="free_y")+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	labs(x="", y="# obs (cells or spots)", title="Pseudobulk sample size")+
	theme_bw()+theme(axis.text.x= element_text(size=7), axis.text.y=element_text(size=8),
		legend.text= element_text(size=7), legend.title=element_text(size=8))

text.df$y_pos = c(4,4)
p2 <- ggplot(full.df, aes(x=x_lab, y=lg10.sum, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=cpList$dx.pal)+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res))+ylim(3,8)+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	labs(x="", y="log10(library size)", title="Library size (# of reads) per pseudobulk sample")+
	theme_bw()+theme(axis.text.x= element_text(size=7))

text.df$y_pos = c(8000,8000)
p3 <- ggplot(full.df, aes(x=x_lab, y=detected, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=cpList$dx.pal)+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res))+ylim(5000,25000)+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	labs(x="", y="# unique genes", title="Detected genes per pseudobulk sample")+
	theme_bw()+theme(axis.text.x= element_text(size=7))

# include prop of transcripts RPS | RPL because this captures the difference between extranuclear and nuclear only
text.df$y_pos = c(10,10)
p4 <- ggplot(full.df, aes(x=x_lab, y=subsets_ribo_percent, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=cpList$dx.pal)+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res))+ylim(0,15)+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	labs(x="", y="percent transcripts RPS | RPL", title="Ribosomal ratio per pseudobulk sample")+
	theme_bw()+theme(axis.text.x= element_text(size=7))



# quantify proportional comparison of number of observations
## make sure to exclude duplicated L2.3 entries so that not inflating the SRT totals
prop.df = filter(full.df, !is.na(clus_orig)) %>% 
#	#aggregate nobs by hand to combine channels
#	group_by(source, condition, sex, clus_orig, sample_id, x_lab, seurat_label, seurat_low.res) %>% summarise(nobs= sum(nobs)) %>%
	# aggregate by dx-sex group
	group_by(source, condition, sex, clus_orig) %>% mutate(group_obs = sum(nobs, na.rm=T)) %>% ungroup() %>%
	mutate(prop_group = nobs/group_obs, source= factor(source, levels=c("MBv","Maitra","SZBD")),
		cond_sex= factor(paste(condition, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")))


#plot annotations
low.res.pal = c("Ast"="#cfa45c","M.V"="#911223",
	"Inb"="#9377AC",
        "L2"="#5D9940", "L3"="#5095CD", 
	"L2.3"=cpList$transfer.bright[["L2.3"]],
        "L4"="#85A0A0",
        "L5"="#ddc94e","L6"="#E45C5F",
        "Olg"="#D1C4B0")

pp1 <- ggplot(distinct(prop.df, x_lab, cond_sex, source, clus_orig, group_obs),
	aes(x=cond_sex, y=group_obs, fill=clus_orig))+
	geom_bar(stat="identity", position="stack", color="grey30", linewidth=.3)+
	scale_fill_manual(values=low.res.pal)+
	facet_grid(cols=vars(source), scales="free_x", space="free_x")+
	labs(y="total # obs (spots or cells)", x="", title="Cell type representation across dx-sex groups")+
	theme_bw()+theme(legend.position="bottom", legend.key.size= unit(8,"pt"))

pp2 <- ggplot(distinct(prop.df, x_lab, cond_sex, source, clus_orig, group_obs),
	aes(x=cond_sex, y=group_obs, fill=clus_orig))+
	geom_bar(stat="identity", position="fill", color="grey30", linewidth=.3)+
	scale_fill_manual(values=low.res.pal)+
	facet_grid(cols=vars(source), scales="free_x", space="free_x")+
	labs(y="prop. obs (spots or cells)", x="", title="How consistent is cell type composition between dx-sex groups?")+
	theme_bw()+theme(legend.position="none")

gridExtra::grid.arrange(pp1, pp2, layout_matrix= matrix(c(1,1,1,2,2), ncol=1))




# calculate per-donor metrics
prop.df = group_by(prop.df, sample_id, x_lab) %>% mutate(donor_obs = sum(nobs)) %>% ungroup() %>%
	mutate(prop_donor= nobs/donor_obs)
# add back in L2.3 duplicate to keep using seurat_low.res facets
prop.df = bind_rows(prop.df, 
	filter(prop.df, source=="MBv", seurat_label=="L2.3") %>% mutate(seurat_low.res="L3", clus_orig=NA)) %>%
	mutate(seurat_low.res = factor(seurat_low.res, levels=levels(sce_mdd$seurat_low.res)))


text.df$y_pos = c(.5,.5)
pp3 <- ggplot(prop.df, aes(x=x_lab, y=prop_donor, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=cpList$dx.pal)+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res))+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	labs(x="", y="prop. obs (spots or cells) per donor", title="How consistent is cell type composition between donors?")+
	theme_bw()+theme(axis.text.x= element_text(size=7))


# hypothetical equal contribution per donor for each dataset depends on the number of donors
# with 20 donors per dataset, equal hypothetical contribution per donor = 1/20 aka .05
# with 10 donors per dataset, equal hypothetical contribution per donor = 1/10 aka 0.1

hypo.df = distinct(full.df, sample_id, x_lab, source, condition, sex) %>% group_by(x_lab, source, condition, sex) %>% tally() %>%
	mutate(equal_prop = 1/n)

text.df$y_pos = c(.4,.4)
pp4 <- ggplot(prop.df, aes(x=x_lab, y=prop_group, color=condition))+
	ggbeeswarm::geom_beeswarm(aes(shape=sex), dodge.width=.8)+
	scale_color_manual(values=c(cpList$dx.pal, "F"="black", "M"="black"))+
	scale_shape_manual(values=c(19,1))+
	facet_wrap(vars(seurat_low.res))+ylim(0,.5)+
	geom_text(data=text.df, aes(x=x_pos, y=y_pos, label=tlabel), color="black", hjust=0)+
	geom_crossbar(data=hypo.df, aes(y=equal_prop, ymin=equal_prop, ymax=equal_prop, color=sex), position = position_dodge(width=.8), width=.5)+
	labs(x="", y="prop. of obs (spots or cells) per dx-sex group", title="How much does any single donor contribution to cell type? (per dx-sex group)",
		subtitle="Hypothetical equal contribution from each donor shown by black bar")+
	theme_bw()+theme(axis.text.x= element_text(size=7))


# save plots
pdf(file="plots/08_compare_DEGs/snRNAseq_QC-abundance_across-studies.pdf", height=7, width=8)
p1
p2
p3
p4
gridExtra::grid.arrange(pp1, pp2, layout_matrix= matrix(c(1,1,1,2,2), ncol=1))
pp3
pp4
dev.off()
cat("\nSaved QC and abundance plots to: plots/08_compare_DEGs/snRNAseq_QC-abundance_across-studies.pdf\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
