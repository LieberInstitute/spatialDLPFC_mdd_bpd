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
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
cat("Dim spe:",dim(spe),"\n")
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "B_1", spe$array)
spe$sample_id2 = paste(spe$slide2, spe$array2)

#subset to batch effect genes plus MOBP for corr and plotting
spe_exclude = spe[c(unique(unlist(exclude.genes)),"ENSG00000168314"),]
cat("Dim spe (batch effect genes + MOBP only):",dim(spe_exclude),"\n")

##in order to reduce the liklihood we are excluding genes that are due to between-sample differences
##in tissue composition (WM biggest concern), correlate all pairs of batch effect genes with MOBP
#### when used with block=spe_exclude, I get an error: "In sqrt((1 - r^2)/(n - 2)) : NaNs produced"
#### probably because most of these are actually batch effect genes and so some values are none for some blocks
#### looked at dotplot.df where prop.expr==0 and mean.expr==0, and the offending genes were sex based and AVP, HSPA6, OXT
#corr.genes = filter(avg.expr.exclude, !gene_name %in% c("XIST","USP9Y","RPS4Y1"))$gene_id # removing just sex-based genes worked
#corr_test = correlatePairs(spe[c(unique(unlist(exclude.genes)),"ENSG00000168314"),], block=spe_exclude$sample_id, equiweight=T, assay.type="logcounts",
#                           subset.row=c(corr.genes,"ENSG00000168314"))
#tmp1 = avg.expr[,c("gene_id","gene_name","avg_expr","decile")]
#tmp2 = avg.expr[,c("gene_id","gene_name","avg_expr","decile")]
#colnames(tmp1) = c("gene1","gene1_name","gene1_avg_expr","gene1_decile")
#colnames(tmp2) = c("gene2","gene2_name","gene2_avg_expr","gene2_decile")
#corr_test = left_join(as.data.frame(corr_test), tmp1) %>% left_join(tmp2)
#write.csv(corr_test, "processed-data/04_feature_selection/batch-effect-genes_with-seq_correlation.csv", row.names=F)
#filter(corr_test, rho>.1)
##few intercorrelations
#filter(corr_test, rho>.1, gene1_name=="MOBP" | gene2_name=="MOBP")

#summarise gene expr for dotplot
spe_summ = scuttle::aggregateAcrossCells(spe_exclude, ids=spe_exclude$sample_id2, 
	statistics=c("mean","prop.detected"),
	use.assay.type="logcounts")


dotplot.df = left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.mean")), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr"),
		tibble::rownames_to_column(as.data.frame(t(scale(t(assay(spe_summ, "logcounts.mean"))))), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="mean_expr_scaled")) %>%
	left_join(tibble::rownames_to_column(as.data.frame(assay(spe_summ, "logcounts.prop.detected")), var="gene_id") %>%
			tidyr::pivot_longer(colnames(spe_summ), names_to="sample_id2", values_to="prop_spots")) %>%
	left_join(avg.expr.exclude)
dotplot.df$y_order = factor(dotplot.df$gene_name, levels=avg.expr.exclude$gene_name[order(avg.expr.exclude$avg_expr)])

#make empty genes for y axis gap
gap.df = filter(dotplot.df, gene_name=="MOBP") %>% mutate(gene_id="gap", mean_expr=0, mean_expr_scaled=0, 
	prop_spots=NA, gene_name="gap", avg_expr=0)
dotplot.df2 = bind_rows(dotplot.df, mutate(gap.df, gene_name="gap1"), mutate(gap.df, gene_name="gap2"), mutate(gap.df, gene_name="gap3"),
	mutate(gap.df, gene_name="gap4"), mutate(gap.df, gene_name="gap5"),mutate(gap.df, gene_name="gap6"), mutate(gap.df, gene_name="gap7"))

#order genes for plotting
all.ordered = c("C3","DDIT4","C5orf63","MTRNR2L1","gap7",
	#iegs
	"NPAS4", "ARC", "NR4A1", "FOS", "DUSP1", "JUNB", "EGR1","gap6",
	#angio-forward stroke group
	"ANGPTL4","VEGFA","CHI3L1","SERPINA3",
	#metallothionein stroke group
	"MT1X",
	#pure stroke genes
	"HAMP","CCL2","ZFP36","C11orf96","GADD45B",
	#oxidative stress stroke group
	"HSPA1B","DNAJB1","HSPA1A","HSPA6","gap5",
	#slide batch
	"AVP","OXT","PURA","PLCG2","ALDOA","gap4",
	#WM/ tissue composition genes (not removing)
	"AQP1", "CERCAM", "NKX6-2","MYRF", "MOBP","gap3",
	#ubiquitous and messy
	"MTRNR2L8","LINC00632","AL627171.2","MALAT1","MAP1B","MTRNR2L12","gap2",
	#L6 and messy
	"PCSK1N","TMSB10","MT3","RPL17","gap1",
	#sex
	"XIST","RPS4Y1","USP9Y")

p1 <- ggplot(dotplot.df2 %>% mutate(y_order3=factor(gene_name, levels=all.ordered)), 
		aes(x=sample_id2, y=y_order3, color=mean_expr_scaled, size=prop_spots))+
	geom_count()+scale_color_viridis_c(option="F", direction=-1)+
	scale_size(range=c(1,4), limits=c(0,1))+labs(color="Avg. expr.\n(scaled)", y="")+
	theme_minimal()+theme(axis.text.x=element_blank())
ggsave("plots/04_feature_selection/batch-effect-genes_dotplot.png", p1, bg="white", height=12, width=16, units="in")
cat("\nSave dotplot to: plots/04_feature_selection/batch-effect-genes_dotplot.png\n")

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

#batch variable plots
cdata2 = group_by(as.data.frame(cdata[,c("sample_id2","slide2","seq","sex","condition","problem_area_flag","sum_umi")]),
                  sample_id2, slide2, seq, sex, condition) %>% summarise(n_flag_spots=sum(problem_area_flag), lg10.total_umi=log10(sum(sum_umi)))

p2 <- ggplot(cdata2, aes(x=sample_id2, y="",fill=n_flag_spots))+
	geom_tile()+scale_fill_gradient(low="#ffff99",high="#b15928")+
	theme_minimal()+theme(axis.text.x=element_blank())
p3 <- ggplot(cdata2, aes(x=sample_id2, y="",fill=lg10.total_umi))+
	geom_tile()+scale_fill_gradient(low="white",high="black")+
	theme_minimal()+theme(axis.text.x=element_blank())

odd_slides = sort(unique(cdata2$slide2))[seq(1,length(unique(cdata2$slide2)), by=2)]
cdata2$slide_bin = ifelse(cdata2$slide2 %in% odd_slides, "odd", "even")
cdata3 = tidyr::pivot_longer(cdata2, c("slide_bin","seq","sex","condition"),
	names_to="coldata", values_to="group")
colpal = c("white","grey","grey20","#619CFF","#F8766D",
	"#b2df8a","#33a02c","#cab2d6","#6a3d9a")
names(colpal) = c("odd","even","NTC","MDD","BPD",
	"SC-TC","SKCCC","Psomagen-1","Psomagen-2")

p4 <- ggplot(filter(cdata3, coldata!="sex") %>% mutate(coldata=factor(coldata, levels=c("seq","condition","slide_bin"))), 
		aes(x=sample_id2, y=coldata, fill=group))+
	geom_tile()+scale_fill_manual(values=colpal)+
	theme_minimal()+theme(axis.text.x=element_blank())

ggsave("plots/04_feature_selection/batch-effect-variables_blocks.png", gridExtra::grid.arrange(p2, p3, p4, ncol=1),
	bg="white", width=12, height=7, units="in")
cat("\nSave heat/tile to: plots/04_feature_selection/batch-effect-variables_blocks.png\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

