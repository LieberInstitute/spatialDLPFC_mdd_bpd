setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scran)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

#load in necessary info
cat("\nNumber of batch effect genes:\n")
exclude.genes = readRDS("processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-seq-sex-condition_list.rds")
length(unique(unlist(exclude.genes))) #47

avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1) %>%
  tibble::rownames_to_column(var="gene_id")
avg.expr.exclude =  filter(avg.expr, gene_id %in% unlist(exclude.genes) | gene_name=="MOBP")

#load in spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
cat("Dim spe:",dim(spe),"\n")
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "B_1", spe$array)
spe$sample_id2 = paste(spe$slide2, spe$array2)

#subset to batch effect genes plus MOBP for corr and plotting
spe_exclude = spe[c(unique(unlist(exclude.genes)),"ENSG00000168314"),]
cat("Dim spe (batch effect genes + MOBP only):",dim(spe_exclude),"\n")

#in order to reduce the liklihood we are excluding genes that are due to between-sample differences
#in tissue composition (WM biggest concern), correlate all pairs of batch effect genes with MOBP
cat("\nCorrelate all batch effect genes	with MOBP\n")
format(Sys.time(), tz="EST")
### when used with block=spe_exclude, I get an error: "In sqrt((1 - r^2)/(n - 2)) : NaNs produced"
### probably because most of these are actually batch effect genes and so some values are none for some blocks
### looked at dotplot.df where prop.expr==0 and mean.expr==0, and the offending genes were sex based and AVP, HSPA6, OXT
corr.genes = filter(avg.expr.exclude, !gene_name %in% c("XIST","USP9Y","RPS4Y1"))$gene_id # removing just sex-based genes worked
corr_test = correlatePairs(spe[c(unique(unlist(exclude.genes)),"ENSG00000168314"),], block=spe_exclude$sample_id, equiweight=T, assay.type="logcounts",
                           subset.row=c(corr.genes,"ENSG00000168314"))
tmp1 = avg.expr[,c("gene_id","gene_name","avg_expr","decile")]
tmp2 = avg.expr[,c("gene_id","gene_name","avg_expr","decile")]
colnames(tmp1) = c("gene1","gene1_name","gene1_avg_expr","gene1_decile")
colnames(tmp2) = c("gene2","gene2_name","gene2_avg_expr","gene2_decile")
corr_test = left_join(as.data.frame(corr_test), tmp1) %>% left_join(tmp2)
write.csv(corr_test, "processed-data/04_feature_selection/batch-effect-genes_correlation.csv", row.names=F)
cat("\nCorrelation results saved to: processed-data/04_feature_selection/batch-effect-genes_correlation.csv\n")
cat("\nCorrelation with MOBP: rho>.1\n")
filter(corr_test, rho>.1, gene1_name=="MOBP" | gene2_name=="MOBP")

#summarise gene expr for dotplot
cat("\n\nSummarise gene expression and perc. expressed per sample\n")
format(Sys.time(), tz="EST")
spe_summ = scuttle::aggregateAcrossCells(spe_exclude, ids=spe_exclude$sample_id2, 
	statistics=c("mean","prop.detected"),
	use.assay.type="logcounts")

cat("\nDotplot\n")
format(Sys.time(), tz="EST")
dotplot.df = left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.mean")), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr"),
		tibble::rownames_to_column(as.data.frame(t(scale(t(assay(spe_summ, "logcounts.mean"))))), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr_scaled")) %>%
	left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.prop.detected")), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="prop_spots")) %>%
	left_join(avg.expr.exclude)
#dotplot.df$y_order = factor(dotplot.df$gene_name, levels=avg.expr.exclude$gene_name[order(avg.expr.exclude$avg_expr)])

#make empty genes for y axis gap
gap.df = filter(dotplot.df, gene_name=="MOBP") %>% mutate(gene_id="gap", mean_expr=0, mean_expr_scaled=0, 
	prop_spots=NA, gene_name="gap", avg_expr=0)
needed.gaps = c(paste0("gap",0:8,".1"),paste0("gap",0:8,".2"))
gap.df = do.call(rbind, lapply(needed.gaps, function(x) mutate(gap.df, gene_name=x)))
dotplot.df2 = bind_rows(dotplot.df, gap.df)

#order genes for plotting
all.ordered = c(#WM/ tissue composition genes (not removing),
  "AQP1", "CERCAM", "NKX6-2","MYRF", "MOBP","gap8.2","gap8.1",
  #iegs
  "NPAS4", "ARC", "NR4A1", "FOS", "DUSP1", "JUNB", "EGR1","gap7.2","gap7.1",
  #angio-forward stroke group
  "ANGPTL4","VEGFA","CHI3L1","SERPINA3","MT1X","gap6.2","gap6.1",
  #pure stroke genes
  "HAMP","CCL2","ZFP36","C11orf96","GADD45B","gap5.2","gap5.1",
  #oxidative stress stroke group
  "HSPA1B","DNAJB1","HSPA1A","HSPA6","gap4.2","gap4.1",
  #sample patterned batch
  "C3","DDIT4","C5orf63","gap3.2","gap3.1",
  #slide patterned
  "AVP","OXT","PURA","PLCG2","ALDOA","gap2.2","gap2.1",
  #decile==10 low spcov
  "MTRNR2L1","PCSK1N","TMSB10","MT3","RPL17","MTRNR2L8","MALAT1","MAP1B","MTRNR2L12",
  #decile==10 high spcov
  "LINC00632","AL627171.2","gap1.2","gap1.1",
  #sex
  "XIST","RPS4Y1","USP9Y","gap0.2","gap0.1")

#add asterisk to indicate if genes qual as SVG and if they would've been removed based on low spcov
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
dotplot.df2$gene_name1 = ifelse(dotplot.df2$gene_name %in% geneList$qual_genes, paste0("*", dotplot.df2$gene_name), dotplot.df2$gene_name)
all.ordered1 = ifelse(all.ordered %in% geneList$qual_genes, paste0("*", all.ordered), all.ordered)
dotplot.df2$gene_name2 = ifelse(dotplot.df2$gene_name %in% setdiff(geneList$qual_genes, geneList$top.decile_low.spcov), 
	paste0("*",dotplot.df2$gene_name1), dotplot.df2$gene_name1)
all.ordered2 = ifelse(all.ordered %in% setdiff(geneList$qual_genes,	geneList$top.decile_low.spcov),	paste0("*",all.ordered1), all.ordered1)

p1 <- ggplot(dotplot.df2 %>% mutate(y_order=factor(gene_name2, levels=all.ordered2)), 
		aes(x=sample_id2, y=y_order, color=mean_expr_scaled, size=prop_spots))+
	geom_count()+scale_color_viridis_c(option="F", direction=-1)+
	scale_size(range=c(1,4), limits=c(0,1))+
	labs(color="Avg. expr.\n(scaled)", y="", subtitle="* = candidate SVG; ** = otherwise considered final SVG")+
	theme_minimal()+theme(axis.text.x=element_blank())
ggsave("plots/04_feature_selection/batch-effect-genes_dotplot.png", p1, bg="white", 
	height=10, width=16, units="in")
cat("\nSave dotplot to: plots/04_feature_selection/batch-effect-genes_dotplot.png\n")

#scatter plot of continuous exp vars
sample.libsize = group_by(as.data.frame(colData(spe)), sample_id) %>% summarise(sum_libsize=sum(sum_umi), med_libsize=median(sum_umi))
summ_cdata = left_join(filter(dotplot.df, gene_name %in% c("MTRNR2L8","LINC00632","AL627171.2","MALAT1","MAP1B","MTRNR2L12","PCSK1N","TMSB10","MT3","RPL17")),
                       as.data.frame(colData(spe))[,c("sample_id2", "sample_id","age", "PMI", "RIN","problem_area_flag")] %>%
                         group_by(sample_id2, sample_id, age, PMI, RIN) %>% summarise(n_flagged_spots=sum(problem_area_flag))) %>% 
  left_join(sample.libsize) %>%
  tidyr::pivot_longer(c("age","PMI","RIN","med_libsize","n_flagged_spots"), names_to="tech_var", values_to="metric") %>%
  mutate(tech_var=factor(tech_var, levels=c("age","PMI","RIN","med_libsize","n_flagged_spots")),
	 gene_name=factor(gene_name, levels=c("RPL17","MT3","TMSB10","PCSK1N","MTRNR2L12","MAP1B","MALAT1","AL627171.2","LINC00632","MTRNR2L8")))

p1 <- ggplot(summ_cdata, aes(x=metric, y=mean_expr_scaled))+
  geom_point(size=.5)+facet_grid(rows=vars(gene_name), cols=vars(tech_var), scales="free_x")+
  theme_minimal()+theme(panel.background = element_rect(color="grey50"))
ggsave("plots/04_feature_selection/batch-effect-genes_continuous-variables_scatter.png",
       p1, bg="white", width=8.5, height=11, units="in")
cat("\nSave dotplot to: plots/04_feature_selection/batch-effect-genes_continuous-variables_scatter.png\n")

#pull in MBv sample info for metadata section of plots (to be combined in AI)
source("code/02_build_spe/getMBvSampleInfo_function.r")
demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")
cdata = merge(colData(spe), demo)
length(unique(cdata$sample_id)) #119
#fix "Mbv_034" to "MBv_034"
cdata[grep("b", cdata$MBv_sample),"MBv_sample"] = "MBv_034"
#fix MBv samples with _DO-NO_SEQ
v1 = unique(cdata$MBv_sample)
names(v1) = v1
v2 = sapply(strsplit(v1,"_"), length)
v2[v2>2]
substr(names(v2[v2>2]), start=0, stop=7)
cdata$MBv_sample = substr(cdata$MBv_sample, start=0, stop=7)
#Samples MBv_001-008 and MBv_013-016 were sequenced on an S4 NovaSeq 6000 (Illumina) at the SC-TC. 
#MBv_017-024 were run on a S4 NovaSeq 6000 at the SKCCC. 
#The remaining samples (MBv_009-012 and MBv_025-120) were sequencing at Psomagen over 9 lanes of a 25B Novaseq X. 
#MBv_121-128, and were sequenced on one lane of a 25B Novaseq X at Psomagen in a fourth sequencing batch
cdata$seq = ifelse(cdata$MBv_sample %in% c(paste0("MBv_00",1:8), paste0("MBv_0",13:16)), "SC-TC", "other")
cdata$seq = ifelse(cdata$MBv_sample %in% paste0("MBv_0",17:24), "SKCCC", cdata$seq)
cdata$seq = ifelse(cdata$MBv_sample %in% c("MBv_009", paste0("MBv_0",10:12), paste0("MBv_0",25:99), paste0("MBv_",100:120)), "Psomagen-1", cdata$seq)
cdata$seq = ifelse(cdata$MBv_sample %in% paste0("MBv_",121:128), "Psomagen-2", cdata$seq)
table(cdata[,c("round","seq")])

#batch variable encoded plots
cdata2 = distinct(as.data.frame(cdata[,c("sample_id2","slide2","seq","sex","condition")]))
odd_slides = sort(unique(cdata2$slide2))[seq(1,length(unique(cdata2$slide2)), by=2)]
cdata2$slide_bin = ifelse(cdata2$slide2 %in% odd_slides, "odd", "even")
cdata3 = tidyr::pivot_longer(cdata2, c("slide_bin","seq","sex","condition"),
	names_to="coldata", values_to="group")
colpal = c("white","grey","grey20","#619CFF","#F8766D",
	"#b2df8a","#33a02c","#cab2d6","#6a3d9a")
names(colpal) = c("odd","even","NTC","MDD","BPD",
	"SC-TC","SKCCC","Psomagen-1","Psomagen-2")

p2 <- ggplot(filter(cdata3, coldata!="sex") %>% mutate(coldata=factor(coldata, levels=c("seq","condition","slide_bin"))), 
		aes(x=sample_id2, y=coldata, fill=group))+
	geom_tile()+scale_fill_manual(values=colpal)+
	theme_minimal()+theme(axis.text.x=element_blank())

ggsave("plots/04_feature_selection/batch-effect-variables_blocks.png", p2,
	bg="white", width=12, height=2.5, units="in")
cat("\nSave batch variable heat/tile to: plots/04_feature_selection/batch-effect-variables_blocks.png\n")

#flag samples for death check
cat("\n\nFlag samples\n")
ieg=c("NPAS4", "ARC", "NR4A1", "FOS", "DUSP1", "JUNB", "EGR1")
bbb.stroke=c("HAMP","CCL2","ZFP36","C11orf96","GADD45B")
angio = c("ANGPTL4","VEGFA","CHI3L1","SERPINA3")
ox.stress = c("HSPA1B","DNAJB1","HSPA1A","HSPA6")

quant.cutoff = filter(dotplot.df, gene_name %in% c(ieg, bbb.stroke, angio, ox.stress)) %>%
  mutate(gene_set=factor(gene_name, levels=c(ieg, ox.stress, angio, bbb.stroke),
                         labels=c(rep("ieg", length(ieg)),
                                  rep("ox.stress", length(ox.stress)),
                                  rep("angio", length(angio)),
                                  rep("bbb.stroke", length(bbb.stroke)))
  )) %>%
  group_by(sample_id2, gene_set) %>% summarise(samp_mean_expr_scaled=mean(mean_expr_scaled)) %>%
  mutate(flag= ifelse(samp_mean_expr_scaled>1, "borderline", "none"),
         flag=ifelse(samp_mean_expr_scaled>2, "flag", flag))

p3 <- ggplot(quant.cutoff, aes(x=sample_id2, y=factor(gene_set, levels=c("ieg","angio","bbb.stroke","ox.stress")), 
                         fill=factor(flag, levels=c("none","borderline","flag"))))+
  geom_tile()+scale_fill_manual("", values=c("white","grey50","black"))+
  labs(title="Per-sample avg. scaled expr. per gene set: >1 = borderline, >2 = flag", y="")+
  theme_minimal()+theme(axis.text.x=element_blank())

ggsave("plots/04_feature_selection/batch-effect_sample-flag_stroke-gene-iegs_blocks.png", p3,
        bg="white", width=12, height=2.5, units="in")
cat("\nSave batch variable heat/tile to: plots/04_feature_selection/batch-effect_sample-flag_stroke-genes-iegs_blocks.png\n")

death.check = filter(quant.cutoff[,c("sample_id2","gene_set","flag")], flag!="none") %>% 
  tidyr::pivot_wider(names_from="gene_set", values_from="flag") %>%
  left_join(distinct(as.data.frame(cdata[,c("sample_id2","brnum","MBv_sample","sample_id")])))
death.check = arrange(death.check[,c("brnum","MBv_sample","sample_id","ieg","ox.stress","angio","bbb.stroke")], MBv_sample)
print(death.check, n=24)
write.csv(death.check, "processed-data/04_feature_selection/check_CoD_stroke-genes-iegs.csv")
cat("\nSave cause of death check table to: processed-data/04_feature_selection/check_CoD_stroke-genes-iegs.csv\n")

#not run
##flag gene set dotplots of logcount expression (not scaled)
##ieg
#ggplot(filter(dotplot.df, gene_name %in% ieg) %>% 
#         mutate(y_order3=factor(gene_name, levels=ieg)),
#       aes(x=sample_id2, y=y_order3, color=mean_expr, size=prop_spots))+
#  geom_count()+scale_color_viridis_c(option="F", direction=-1)+#, limits=c(0,1.8))+
#  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.", y="")+
#  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
#                        legend.position="bottom")

##bbb.stroke
#ggplot(filter(dotplot.df, gene_name %in% bbb.stroke) %>% 
#        mutate(y_order3=factor(gene_name, levels=bbb.stroke)),
#      aes(x=sample_id2, y=y_order3, color=mean_expr, size=prop_spots))+
#  geom_count()+scale_color_viridis_c(option="F", direction=-1, limits=c(0,1))+
#  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.", y="")+
#  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
#                        legend.position="bottom")

##ox.stress
#ggplot(filter(dotplot.df, gene_name %in% ox.stress) %>% 
#        mutate(y_order3=factor(gene_name, levels=ox.stress)),
#      aes(x=sample_id2, y=y_order3, color=mean_expr, size=prop_spots))+
#  geom_count()+scale_color_viridis_c(option="F", direction=-1)+#, limits=c(0,1.8))+
#  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.", y="")+
#  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
#                        legend.position="bottom")

##angio
#ggplot(filter(dotplot.df, gene_name %in% angio) %>% 
#         mutate(y_order3=factor(gene_name, levels=angio)),
#       aes(x=sample_id2, y=y_order3, color=mean_expr, size=prop_spots))+
#  geom_count()+scale_color_viridis_c(option="F", direction=-1)+#, limits=c(0,1.8))+
#  scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.", y="")+
#  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
#                        legend.position="bottom")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

