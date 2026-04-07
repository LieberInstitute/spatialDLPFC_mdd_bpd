suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
  library(gridExtra)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
cat("\nHow many spots are dropped after smoothing clusters:\n")
table(cdata[cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),"smoothed_k9_1663_f"])

cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]

spe = spe[,rownames(cdata2)]
spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM"),
                              labels=c("L1","L2","L3.4","L5","L6","WM"))

#load aucell for checking inhb (b/c pop of interest) and astro (b/c large sample-to-sample variability)
aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)],
                                                                        function(x) substr(x, start=0, stop=nchar(x)-3))

seurat_levels= c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")
cond_sex = c("F NTC","F MDD", "F BPD","M NTC","M MDD","M BPD")

aucell$sample_id = substr(rownames(aucell), start=20, stop=50)
aucell$sample_id2 = paste(aucell$sample_id, aucell$seurat_label)
aucell$seurat_label= factor(aucell$seurat_label, levels=c("Micro.Vasc", "Astro", "L2.3", "L4", "Inhb", "L5", "L6", "Oligo"),
                            labels=seurat_levels)
aucell$sex.group = paste(aucell$sex, aucell$condition)
aucell$clus_groups= factor(paste(aucell$sex.group, aucell$seurat_label),
                           levels=as.character(outer(cond_sex, seurat_levels, paste)))


#all example samples that are NTC
rnascope.samples = c("342-A1", "332-A1", "279-A1", "327-C1", #NTC F
                     "308-D1", "329-A1", "342-B1", "332-B1") #NTC M
rnascope.samples = paste0("V13B23-", gsub("-","_", rnascope.samples))

dim(aucell)
colnames(aucell)
aucell2 = filter(aucell, sample_id %in% rnascope.samples)

col.pal = cpList$transfer.bright
names(col.pal) = c("M.V","Ast","L2.3","L4","L5","L6","Olg","Inb")
p1 <- ggplot(aucell2, aes(x=sample_id, fill=seurat_label))+
  geom_bar(stat="count", position="fill")+facet_wrap(vars(sex.group), scales="free_x")+
  scale_fill_manual("",values=col.pal)+theme_minimal()+
  theme(axis.text.x=element_text(size=8, angle=45, hjust=1))
p2 <- ggplot(aucell2, aes(x=sample_id, fill=smoothed_k9_1663))+
  geom_bar(stat="count", position="fill")+facet_wrap(vars(sex.group), scales="free_x")+
  scale_fill_manual("", values=cpList$smoothed.bright)+theme_minimal()+
  theme(axis.text.x=element_text(size=8, angle=45, hjust=1))

grid.arrange(p1, p2, ncol=1)


p1 <- ggplot(aucell2, aes(x=LHX6, group=sample_id2, color=seurat_label))+
  stat_ecdf()+scale_color_manual(values=col.pal, guide="none")+theme_minimal()+theme(aspect.ratio=1)

p2 <- ggplot(aucell2, aes(x=ARX, group=sample_id2, color=seurat_label))+
  stat_ecdf()+scale_color_manual(values=col.pal, guide="none")+theme_minimal()+theme(aspect.ratio=1)

#p3 <- ggplot(aucell2, aes(x=DLX1, group=sample_id2, color=seurat_label))+
#  stat_ecdf()+scale_color_manual(values=col.pal, guide="none")+theme_minimal()+theme(aspect.ratio=1)
p3 <- ggplot(aucell2, aes(x=SOX2, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

p4 <- ggplot(aucell2, aes(x=SOX9, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

grid.arrange(p1, p2, p3, p4, ncol=2)


filter(aucell2, seurat_label=="Inb") %>% group_by(sex.group, sample_id) %>% 
         summarise(med_LHX6=median(LHX6))
group_by(aucell2, sex.group, sample_id) %>% 
  summarise(med_SOX9=median(SOX9))



cdata3 = filter(as.data.frame(colData(spe)), sample_id %in% rnascope.samples)

p1 <- ggplot(filter(aucell2, seurat_label=="Inb"), 
             aes(x=LHX6, group=sample_id, lty=sex))+
  stat_ecdf()+labs(x="LHX6 (Inhb spots only)")+theme_minimal()+theme(aspect.ratio=1)

p2 <- ggplot(aucell2, aes(x=SOX9, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

p3 <- ggplot(cdata3, aes(x=sum_umi, group=sample_id, lty=sex))+
  stat_ecdf()+scale_x_log10()+theme_minimal()+theme(aspect.ratio=1)

p4 <- ggplot(cdata3, aes(x=sum_gene, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

grid.arrange(p1, p2, p3, p4, ncol=2)

group_by(cdata3, condition, sex, sample_id) %>%
  summarise(med_detected=median(sum_gene))


#rnascope.samples = c("342-A1", "279-A1", #NTC F
#                     "329-A1", "342-B1") #NTC M
rnascope.samples = c("342-A1", "332-A1", #NTC F
                     "308-D1", "342-B1") #NTC M

rnascope.samples = paste0("V13B23-", gsub("-","_", rnascope.samples))

p1 <- ggplot(filter(aucell2, sample_id %in% rnascope.samples), aes(x=sample_id, fill=seurat_label))+
  geom_bar(stat="count", position="fill")+facet_wrap(vars(sex.group), scales="free_x")+
  scale_fill_manual("",values=col.pal)+theme_minimal()#+
  #theme(axis.text.x=element_text(size=8, angle=45, hjust=1))
p2 <- ggplot(filter(aucell2, sample_id %in% rnascope.samples), aes(x=sample_id, fill=smoothed_k9_1663))+
  geom_bar(stat="count", position="fill")+facet_wrap(vars(sex.group), scales="free_x")+
  scale_fill_manual("", values=cpList$smoothed.bright)+theme_minimal()#+
  #theme(axis.text.x=element_text(size=8, angle=45, hjust=1))

grid.arrange(p1, p2, ncol=1)

p1 <- ggplot(filter(aucell2, seurat_label=="Inb", sample_id %in% rnascope.samples), 
             aes(x=LHX6, group=sample_id, lty=sex))+
  stat_ecdf()+labs(x="LHX6 (Inhb spots only)")+theme_minimal()+theme(aspect.ratio=1)

p2 <- ggplot(filter(aucell2, sample_id %in% rnascope.samples), aes(x=SOX9, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

p3 <- ggplot(filter(cdata3, sample_id %in% rnascope.samples), aes(x=sum_umi, group=sample_id, lty=sex))+
  stat_ecdf()+scale_x_log10()+theme_minimal()+theme(aspect.ratio=1)

p4 <- ggplot(filter(cdata3, sample_id %in% rnascope.samples), aes(x=sum_gene, group=sample_id, lty=sex))+
  stat_ecdf()+theme_minimal()+theme(aspect.ratio=1)

grid.arrange(p1, p2, p3, p4, ncol=2)


spot.genes = c("GAD1","TAC1","TRBC2",
               "GFAP","AIF1","MTRNR2L1")

spe_tmp = spe[,spe$sample_id %in% rnascope.samples]

for(i in spot.genes) {
  spe_tmp[[i]] = logcounts(spe_tmp)[rowData(spe_tmp)$gene_name==i,]
}


spotList = list()


for(y in spot.genes) {
  max.val = max(colData(spe_tmp)[[y]])
  max.valr = round(max.val,1)
  if(max.valr<max.val) max.valr=max.valr+.1
  
  plist1 <- lapply(rnascope.samples, function(x) {
    spe_sub = spe_tmp[,spe_tmp$sample_id==x]
    p = make_escheR(spe_sub) |> add_ground(var="smoothed_k9_1663", 
                                           stroke=.3, point_size = .5) |> 
      add_fill(var=y, point_size = .5)
    p2 <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
      scale_fill_gradient(limits=c(0,max.valr), low="white",high="black", guide="none")+
      labs(subtitle=unique(spe_sub$sample_id))+
      theme(plot.subtitle=element_text(size=9))
    return(rasterize(p2, dpi=200))
  })
  
  spotList[[y]] = arrangeGrob(grobs=plist1, ncol=4, nrow=1, top=paste(y, "(fill scale fixed, can compare across sections)"))
}

ggsave(file="plots/publication/Figure3/rnascope-samples_key-spot-plots.pdf", 
       marrangeGrob(grobs=spotList, ncol=1, nrow=6, top=NULL),
       height=12, width=9)
