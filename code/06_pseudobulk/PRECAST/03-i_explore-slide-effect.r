# JT: see locally saved file `explore_slide_effect.r` for more plots, exploration, and categorical consideration of new batch variable

library(SpatialExperiment)
library(edgeR)
library(scater)
library(dplyr)
library(ggplot2)

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")

exp.vars = c("precast_k9_1663","sample_id",
             "slide","round",
             "condition","sex",
             "age","PMI","RIN",
             "sum","detected","nspots",
             "subsets_mito_percent", "subsets_ribo_percent")

var.m = getVarianceExplained(spe_pseudo, variables=exp.vars)

#checked 40, 50, and 60
slide.genes_40 = rownames(var.m)[var.m[,"slide"]>40] 
slide.genes_50 = rownames(var.m)[var.m[,"slide"]>50] 
slide.genes_60 = rownames(var.m)[var.m[,"slide"]>60] 

#add vars of experimenter and date for heatmap column annotation
REDCap <- read.csv("raw-data/sample_info/Visium_DATA_2025-01-22_1406.csv")
{
  A1 <- subset(REDCap, select = c("slide", "sample_a1", "project_a1","sample_number1_a1","date","experimenter"))
  B1 <- subset(REDCap, select = c("slide", "sample_b1", "project_b1","sample_number1_b1","date","experimenter"))
  C1 <- subset(REDCap, select = c("slide", "sample_c1", "project_c1","sample_number1_c1","date","experimenter"))
  D1 <- subset(REDCap, select = c("slide", "sample_d1", "project_d1","sample_number1_d1","date","experimenter"))
  colnames(A1) <- colnames(B1) <- colnames(C1) <- colnames(D1) <- c("slide", "sample", "project","MBv_sample","date","experimenter")
  A1$array <- "A1"
  B1$array <- "B1"
  C1$array <- "C1"
  D1$array <- "D1"
  REDCap_table <- rbind(A1, B1, C1, D1)
  REDCap_table <- REDCap_table[order(REDCap_table$slide), ]
  
  #filter to mbv
  REDCap_MBv <- REDCap_table[which(REDCap_table$project == "spatialDLPFC_MBv_4100"), ]
  
  #fix slide 301 array info (see note in raw-data/sample_info/README
  REDCap_MBv[REDCap_MBv[["slide"]]=="V13B23-301","array"] = c("D1","C1","B1","A1")
  
  #fix "Mbv_034" to "MBv_034"
  REDCap_MBv[grep("b", REDCap_MBv[["MBv_sample"]]),"MBv_sample"] = "MBv_034"
  
  #fix MBv samples with _DO-NO_SEQ
  v1 = unique(REDCap_MBv[["MBv_sample"]])
  names(v1) = v1
  v2 = sapply(strsplit(v1,"_"), length)
  #v2[v2>2]
  #substr(names(v2[v2>2]), start=0, stop=7)
  REDCap_MBv[["MBv_sample"]] = substr(REDCap_MBv[["MBv_sample"]], start=0, stop=7)
}

REDCap_MBv[REDCap_MBv[["MBv_sample"]] %in% spe_pseudo$MBv_sample,]
colData(spe_pseudo) = merge(colData(spe_pseudo), REDCap_MBv[REDCap_MBv[["MBv_sample"]] %in% spe_pseudo$MBv_sample,])

#plot heatmap for visualization and for pulling out hclust of rows
tmp = distinct(as.data.frame(colData(spe_pseudo)), sample_id, slide, experimenter, date)
col_annot = data.frame("slide"=tmp$slide, row.names = tmp$sample_id)
col_annot$experimenter = gsub("Stephanie Page", "SCP", tmp$experimenter)
col_annot$experimenter = gsub("Lina Oh", "Lina", col_annot$experimenter)
col_annot$experimenter = gsub("Svitlana Bach", "SB", col_annot$experimenter)
table(col_annot$experimenter)
col_annot$date = tmp$date

col.list = unlist(sapply(1:length(unique(tmp$slide)), function(x) {
  u.slide = unique(tmp$slide)
  v1 = u.slide[[x]]
  if(round(x/2)==(x/2)) {
    col1 = "black"
  } else {col1 = "grey80"}
  names(col1) = v1
  return(col1)
}))
annot_colors=list("slide"=col.list)
annot_colors$experimenter = RColorBrewer::brewer.pal(n=4, "Dark2")
names(annot_colors$experimenter) = unique(col_annot$experimenter)

#colors for date
as.Date(col_annot$date)
table(as.Date(col_annot$date))

annot_colors$date = c("#E2B1C6","#FF8A97","#E872AA","#941D61", #"oct 16","oct 23","nov 1","dec 18",
  "#98DA9C","#7ACB75", #Feb 5,6
  "#599AD7","#3495CA","#2C87A1", #Feb 12,13,15
  "#BCD5ED","#96C6E1","#71BDD5",#Feb 20,21,22
  "#9EB3B6","#BAC8CD","#D6DEE3",#Feb 26,29, Mar 1
  "#D7F4AC","#D1E457","#DBDA2F",#Mar 12,14,15
  "#FAD590","#F5B876","#F09B5D","#EC7F43","#E7622A","#E34611", #Mar 20,22,23,25,26,27
  "#A76E51","#BD936F","#DFD2B1",#Apr 1,3,5
  "#463929",#Apr 18
  "#D7A5EE","#8657BF" #Nov 8, 11
  )
names(annot_colors$date) <- as.character(sort(unique(as.Date(col_annot$date))))

hm_40 = plotGroupedHeatmap(spe_pseudo, features= slide.genes_40, 
                           group = "sample_id", cluster_cols=F, center=T,
                           annotation_col=col_annot, annotation_colors=annot_colors, 
                           show_rownames=F)#, silent=T)

#heatmap in date order
{
  spe_pseudo$sample_factor = factor(paste(spe_pseudo$date, spe_pseudo$sample_id),
                                    levels=sort(paste(col_annot$date, rownames(col_annot))),
                                    labels=substr(sort(paste(col_annot$date, rownames(col_annot))), 12, 24))
  plotGroupedHeatmap(spe_pseudo, #features= slide.genes_60, 
                     features= rowData(spe_pseudo)[slide.genes_60,"gene_name"], swap_rownames="gene_name", 
                     group = "sample_factor", cluster_cols=F, center=T,
                     annotation_col=col_annot, annotation_colors=annot_colors, 
                     #show_rownames=F, 
                     fontsize_col=8)
  plotGroupedHeatmap(spe_pseudo, features= both.slide.genes, swap_rownames="gene_name", 
                     group = "sample_factor", cluster_cols=F, center=T,
                     annotation_col=col_annot, annotation_colors=annot_colors, 
                     fontsize_col=8)
}


#split genes intro groups
gene.groups_40 = cutree(hm_40$tree_row, k=2)
table(gene.groups_40)
#more 40
#. 1   2 
#220 118
gene.group1_40 = names(gene.groups_40)[gene.groups_40==1]
gene.group2_40 = names(gene.groups_40)[gene.groups_40==2]

#aggregate gex by group
cdata = as.data.frame(colData(spe_pseudo))
cdata$more40_gene.group1 = colSums(logcounts(spe_pseudo)[gene.group1_40,])
cdata$more40_gene.group2 = colSums(logcounts(spe_pseudo)[gene.group2_40,])


#extra visuals
ggplot(cdata, aes(x=subsets_ribo_percent, y=more40_gene.group1))+
  geom_point()+facet_wrap(vars(precast_k9_1663))+
  theme_bw()+ggtitle(">40% var. to slide")

ggplot(cdata, aes(x=subsets_ribo_percent, y=more40_gene.group2))+
  geom_point()+facet_wrap(vars(precast_k9_1663))+
  theme_bw()+ggtitle(">40% var. to slide")

#collapse to sample level
s.cdata = group_by(cdata, sample_id) %>% 
  summarise(avg_more40_group1=mean(more40_gene.group1),
            avg_more40_group2=mean(more40_gene.group2),
           # avg_more50_group1=mean(more50_gene.group1),
           # avg_more50_group2=mean(more50_gene.group2),
           # avg_more60_group1=mean(more60_gene.group1),
           # avg_more60_group2=mean(more60_gene.group2),
            age=unique(age), PMI=unique(PMI), RIN=unique(RIN), seq=unique(seq),
            condition=unique(condition), sex=unique(sex), slide=unique(slide))

cor(scale(s.cdata$avg_more40_group1), scale(s.cdata$avg_more40_group2))
#-0.8222957


#generate continuous batch effect variable at the level of capture area
m40 = lm(scale(s.cdata$avg_more40_group2) ~ scale(s.cdata$avg_more40_group1))
summary(m40) 
s.cdata$m40_fitted = m40$fitted.values
ggplot(s.cdata, aes(x=scale(avg_more40_group1), y=scale(avg_more40_group2), color=m40_fitted))+
  geom_point()+geom_abline(aes(intercept=coef(m40)[[1]], slope=coef(m40)[[2]]))+
  theme_bw()+labs(title=">40% var. to slide", x="group 1 genes", y="group2 genes")


#check for conflation with other known variables
#plot against continuous variables
#age, pmi, rin
p1 <- ggplot(s.cdata, aes(x=m40_fitted, y=age))+
  geom_point()+theme_bw()
p2 <- ggplot(s.cdata, aes(x=m40_fitted, y=PMI))+
  geom_point()+theme_bw()
p3 <- ggplot(s.cdata, aes(x=m40_fitted, y=RIN))+
  geom_point()+theme_bw()
#seq, condition, sex
p4 <- ggplot(s.cdata, aes(x=seq, y=m40_fitted))+
  ggbeeswarm::geom_quasirandom()+geom_boxplot(fill="white", alpha=.5)+
  theme_bw()
p5 <- ggplot(s.cdata, aes(x=condition, y=m40_fitted))+
  ggbeeswarm::geom_quasirandom()+geom_boxplot(fill="white", alpha=.5)+
  theme_bw()
p6 <- ggplot(s.cdata, aes(x=sex, y=m40_fitted))+
  ggbeeswarm::geom_quasirandom()+geom_boxplot(fill="white", alpha=.5)+
  theme_bw()

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)

ggplot(s.cdata, aes(x=m40_fitted, color=condition))+
  geom_density(linewidth=2)+#stat_ecdf(linewidth=2)+
  scale_color_manual(values=cpList$dx.pal)+
  facet_wrap(vars(sex))+
  theme_bw()

#slide, 
tmp = group_by(s.cdata, slide) %>% summarise(m40_avg = mean(m40_fitted))
tmp$slide_order = rank(tmp$m40_avg)
s.cdata$slide_order = factor(s.cdata$slide, levels=tmp$slide[order(tmp$slide_order)])
ggplot(s.cdata, aes(x=slide_order, y=m40_fitted, color=condition))+
  geom_point()+geom_boxplot(fill="white", color="black", alpha=.5)+
  scale_color_manual(values=cpList$dx.pal)+theme_bw()

#tissue composition
s.cdata2 = group_by(cdata, sample_id) %>% mutate(total_spots=sum(nspots)) %>% 
  filter(precast_k9_1663=="L1") %>% mutate(perc.L1 = nspots/total_spots) %>%
  select(sample_id, perc.L1) %>% right_join(s.cdata)

s.cdata2 = group_by(cdata, sample_id) %>% mutate(total_spots=sum(nspots)) %>% 
  ungroup() %>% mutate(perc_spots=nspots/total_spots) %>% select(sample_id, precast_k9_1663, perc_spots) %>%
  tidyr::pivot_wider(names_from="precast_k9_1663", names_prefix="perc_", values_from="perc_spots", values_fill=0) %>%
  left_join(s.cdata[,c("sample_id","m40_fitted")]) %>%
  tidyr::pivot_longer(all_of(paste0("perc_",c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))),
                      names_to="cluster", names_prefix = "perc_", values_to="perc_spots")
head(s.cdata2)
ggplot(s.cdata2, aes(x=m40_fitted, y=perc_spots))+
  geom_point()+facet_wrap(vars(cluster), scale="free_y")+
  theme_bw()


#save sample-level
write.csv(s.cdata, "processed-data/06_pseudobulk/continuous_batch_variable_slide-sample-id.csv", row.names=F)

