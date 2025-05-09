#script was run in interactive session
library(SpatialExperiment)
library(ggplot2)
library(dplyr)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
cdata = as.data.frame(colData(spe_pseudo))

source("code/05_clustering/PRECAST/PRECAST_colorLists.r")

#candidate sections for vistoseg
#308, 309, 327, 329, 332, 333 (A1 is weirdo massive L5), 342, 352, 380, 381, 382
filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  group_by(condition, sex) %>% tally()

#great! enough to choose from so pick top 4 from each 
filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="NTC", sex=="M") %>% arrange(slide)
#381-C1, 342-B1, 332-B1, 329-A1, 308-D1
ntc.m = c("381-C1", "342-B1", "332-B1", "329-A1", "308-D1")

filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="NTC", sex=="F") %>% arrange(slide)
#327-C1 (smaller side), 329-B1, 332-A1, 342-A1
#308-C1 missing some L1 but not much
ntc.f = c("327-C1", "329-B1", "332-A1", "342-A1")

filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="MDD", sex=="M") %>% arrange(slide)
#382-C1, 380-A1 (funky but good), 352-A1, 329-C1, 309-D1
#333-B1 (no WM), 
mdd.m = c("382-C1", "380-A1", "352-A1", "329-C1", "309-D1")

filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="MDD", sex=="F") %>% arrange(slide)
#382-D1, 380-B1, 352-B1, 329-D1, 309-C1
mdd.f = c("382-D1", "380-B1", "352-B1", "329-D1", "309-C1")

filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="BPD", sex=="M") %>% arrange(slide)
#382-A1, 352-C1, 342-D1, 332-D1 (funky but good), 327-B1, 
#381-A1 (lots of WM but good), 333-D1 (avoid 333 if possible bc low umi)
bpd.m = c("382-A1", "352-C1", "342-D1", "332-D1", "327-B1")


filter(cdata, slide %in% paste("V13B23",c(308, 309, 327, 329, 332, 333, 342, 352, 380, 381, 382), sep="-")) %>%
  distinct(sample_id, slide, array, condition, sex) %>% 
  filter(condition=="BPD", sex=="F") %>% arrange(slide)
#382-B1, 352-D1 (higher side of Vasc), 332-C1, 327-A1 (small but good), 309-A1, 308-A1
#381-B1 (higher WM but good), 342-C1 (lots of L2 but good), 333-C1 (avoid 333 but good)
bpd.f = c("382-B1", "352-D1", "332-C1", "327-A1", "309-A1", "308-A1")

candidate.samples = paste0("V13B23-", gsub("-","_",c(ntc.m, ntc.f, mdd.m, mdd.f, bpd.m, bpd.f)))

ggplot(filter(cdata, sample_id %in% candidate.samples), aes(x=sample_id, y=ncells, fill=combined_cluster))+
  geom_bar(stat="identity", position="stack", color="black", linewidth=.3)+facet_grid(cols=vars(sex), rows=vars(condition), scales="free_x")+
  scale_fill_manual(values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

ggplot(filter(cdata, sample_id %in% candidate.samples), aes(x=sample_id, y=ncells, fill=combined_cluster))+
  geom_bar(stat="identity", position="fill", color="black", linewidth=.3)+facet_grid(cols=vars(sex), rows=vars(condition), scales="free_x")+
  scale_fill_manual(values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


#final samples
vistoseg.samples = c("342-B1", "332-B1", "329-A1", "308-D1", #NTC M
                     "327-C1", "329-B1", "332-A1", "342-A1", #NTC F
                     "382-C1", "352-A1", "329-C1", "309-D1", #MDD M
                     "382-D1", "352-B1", "309-C1", "380-B1", #"329-D1", #MDD F
                     "382-A1", "352-C1", "342-D1", "327-B1", #BPD M
                     "382-B1", "327-A1", "309-A1", "308-A1") #BPD F

vistoseg.samples = paste0("V13B23-", gsub("-","_", vistoseg.samples))

cdata2 = filter(cdata, sample_id %in% vistoseg.samples) %>% 
  mutate(cond_sex= factor(paste(condition, sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F")))

ggplot(cdata2, aes(x=sample_id, y=ncells, fill= combined_cluster))+
  geom_bar(stat="identity", position="fill", color="black", linewidth=.3)+
  facet_wrap(vars(cond_sex), scales="free_x", ncol=6)+
  scale_fill_manual(values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+theme(axis.text.x= element_text(angle=90, hjust=1, vjust=.5))



#spot plots
library(HDF5Array)
library(ggspavis)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
dim(spotdata) # 535248     19
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))

spotdata = spotdata[!is.na(spotdata$combined_cluster),]
dim(spotdata) # 518935     19
spe = spe[,rownames(spotdata)]
dim(spe) # 28965 518935
spe$combined_cluster = factor(spotdata$combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
table(spe$combined_cluster, useNA="ifany")


spe_sub = spe[,spe$sample_id %in% vistoseg.samples]
dim(spe_sub)
# 28965 110922 #with 380-B1
# 28965 111290 #with 329-D1
spe_sub$cond_sex = factor(paste(spe_sub$condition, spe_sub$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))
spe_sub$facet_col = factor(spe_sub$sample_id, levels=vistoseg.samples, labels=paste0("C",rep(c(1:4), 6)))

mod_spatialCoords = spatialCoords(spe_sub)
for (i in unique(spe_sub$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}
dim(mod_spatialCoords)

#find coordinates where all samples have spot to use for ifelse encoding to place sample_id only once for each facet
group_by(as.data.frame(colData(spe_sub)), array_row, array_col) %>% tally() %>% filter(n==24)
#array_row = 5, array_col = 13

plotSpots(spe_sub, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
          sample_id="sample_id", annotate="combined_cluster", point_size=.1)+
  scale_color_manual("combined\ncluster", values=precast.colorList[["n1663_k9"]]$colors[1:8])+
  facet_grid(rows=vars(cond_sex), cols=vars(facet_col), switch="y")+
  geom_text(aes(x=0, y= -60, 
                label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13,
                              .data[["sample_id"]], "")), hjust=0, vjust=0, 
            color="black", size=2)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

#get sequencing round information for MBv sample number
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

df = filter(as.data.frame(cdata), sample_id %in% vistoseg.samples) %>% 
  select(sample_id, brnum, MBv_sample, condition, sex, age, PMI, RIN) %>%
  distinct()

write.csv(df, "raw-data/sample_info/Vistoseg_samples_n24.csv", row.names=F)

