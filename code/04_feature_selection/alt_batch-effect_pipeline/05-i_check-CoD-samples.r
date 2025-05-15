library(SpatialExperiment)
library(HDF5Array)
library(DelayedArray)
library(scran)
library(scater)
library(dplyr)
library(ggplot2)

set.seed(123)

avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1) %>%
  tibble::rownames_to_column(var="gene_id")
exclude.genes = readRDS("processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-seq-sex-condition_list.rds")
#old.exclude.genes = readRDS("processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-sex-condition_list.rds")
length(unique(unlist(exclude.genes))) #47

avg.expr.exclude =  filter(avg.expr, gene_id %in% unlist(exclude.genes) | gene_name=="MOBP")
#avg.expr.exclude =  filter(avg.expr, gene_id %in% unlist(exclude.genes) | gene_name %in% c("GFAP","MOBP"))

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
#spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "B_1", spe$array)
spe$sample_id2 = paste(spe$slide2, spe$array2)
sample.libsize = group_by(as.data.frame(colData(spe)), sample_id) %>% summarise(sum_libsize=sum(sum_umi), med_libsize=median(sum_umi))

#add MOBP for correlation
#MOBP: "ENSG00000168314"
#GFAP: "ENSG00000131095"
spe_exclude = spe[c(unique(unlist(exclude.genes)),"ENSG00000168314"),]

#dotplot
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


#get sequencing round information
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
length(unique(cdata$sample_id2)) #119

#flag samples
#just plot flag gene lists with true avg expr
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

ggplot(quant.cutoff, aes(x=sample_id2, y=factor(gene_set, levels=c("ieg","angio","bbb.stroke","ox.stress")), 
                         fill=factor(flag, levels=c("none","borderline","flag"))))+
  geom_tile()+scale_fill_manual("", values=c("white","grey50","black"))+
  labs(title="Per-sample avg. scaled expr. per gene set: >1 = borderline, >2 = flag", y="")+
  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
                        legend.position="bottom")

death.check = filter(quant.cutoff[,c("sample_id2","gene_set","flag")], flag!="none") %>% 
  tidyr::pivot_wider(names_from="gene_set", values_from="flag") %>%
  left_join(distinct(as.data.frame(cdata[,c("sample_id2","brnum","MBv_sample","sample_id")])))
death.check = arrange(death.check[,c("brnum","MBv_sample","sample_id","ieg","ox.stress","angio","bbb.stroke")], MBv_sample)
print(death.check, n=24)
write.csv(death.check, "processed-data/04_feature_selection/check_CoD_stroke-genes-iegs.csv")
