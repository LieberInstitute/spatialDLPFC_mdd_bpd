setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
dim(spe_pseudo) 
spe_sm <- spe_pseudo

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
dim(spe_pseudo) 
spe_se <- spe_pseudo


# plot 1: pc3 gene loadings b/w annotation strategies
rot.mtx_sm = attr(reducedDim(spe_sm, "PCA_1663"), "rotation")
rot.mtx_se = attr(reducedDim(spe_se, "PCA_1663"), "rotation")

r1 = data.frame("pc3_sm"=rot.mtx_sm[rownames(rot.mtx_sm),3], 
                "pc3_se"=rot.mtx_se[rownames(rot.mtx_sm),3],
                row.names=rownames(rot.mtx_sm))
r1$gene_name = rowData(spe_sm)[rownames(r1),"gene_name"]

genes1 = c("MT1X","RPS26","HBA1","HBA2",
  "AL627171.2","ZBTB20","PLCG2",
  "MALAT1","MTRNR2L8","MTRNR2L12")

p1 <- ggplot(r1, aes(x=pc3_sm, y=pc3_se))+
  geom_point(size=.3, color="grey")+
  ggrepel::geom_text_repel(data=filter(r1, gene_name %in% genes1),
                            aes(label=gene_name), size=3, min.segment.length=0,
                           fontface="italic")+
  labs(x="PRECAST (smoothed) gene loadings", y="Seurat PC30 gene loadings",
       title="PC3 rotation for n=1663 SVGs")+
  theme_minimal()+theme(text=element_text(size=9))


# plot 2: slide variance explained b/w annotation strategies
spe_sm$pc3 = reducedDim(spe_sm, "PCA_1663")[,3]
spe_se$pc3 = reducedDim(spe_se, "PCA_1663")[,3]

var.m_sm = getVarianceExplained(spe_sm, variables=c("smoothed_k9_1663","sample_id","slide","seq","pc3"))
w.max = apply(var.m_sm, MARGIN=1, which.max)
table(colnames(var.m_sm)[w.max])
#sample_id smoothed_k9_1663 
#9101             4743

var.m_se = getVarianceExplained(spe_se, variables=c("seurat_label","sample_id","slide","seq","pc3"))
w.max = apply(var.m_se, MARGIN=1, which.max)
table(colnames(var.m_se)[w.max])
#pc3    sample_id seurat_label 
#2         8854         4263


r2 = data.frame("slide_sm"=var.m_sm[rownames(rot.mtx_sm),"slide"], 
                "slide_se"=var.m_se[rownames(rot.mtx_sm),"slide"],
                row.names=rownames(rot.mtx_sm))
r2$gene_name = rowData(spe_sm)[rownames(r2),"gene_name"]


p2 <- ggplot(r2, aes(x=slide_sm, y=slide_se))+
  geom_point(size=.3, color="grey")+
  ggrepel::geom_text_repel(data=filter(r2, gene_name %in% genes1),
                           aes(label=gene_name), size=3, min.segment.length=0,
                           fontface="italic")+
  labs(x="PRECAST (smoothed) gene variance", y="Seurat PC30 gene variance",
       title="Variance explained by slide for n=1663 SVGs")+
  theme_minimal()+theme(text=element_text(size=9))


# plots 3 and 4: pc3 gene loadings vs slide variance explained
d1_sm = data.frame("pc3"=rot.mtx_sm[,"PC3"], "var_pc3"=var.m_sm[rownames(rot.mtx_sm),"pc3"],
                "var_slide"=var.m_sm[rownames(rot.mtx_sm),"slide"],
                row.names = rownames(rot.mtx_sm))
d1_sm$gene_name = rowData(spe_sm)[rownames(d1_sm),"gene_name"]

d1_se = data.frame("pc3"=rot.mtx_se[,"PC3"], "var_pc3"=var.m_se[rownames(rot.mtx_se),"pc3"],
                   "var_slide"=var.m_se[rownames(rot.mtx_se),"slide"],
                   row.names = rownames(rot.mtx_se))
d1_se$gene_name = rowData(spe_se)[rownames(d1_se),"gene_name"]

genes2 = c("APC","MTRNR2L10","LINC00632","ANKRD12","HMGB1","NTM","DST",
           "ELOB","GUK1","MRPL41","MT3","NDUFA3","MAP1LC3A","COX4I1")

p3<- ggplot(d1_sm, aes(x=pc3, y=var_slide))+
  geom_point(size=.3, color="grey")+
  ggrepel::geom_text_repel(data=filter(d1_sm, gene_name %in% genes2),
                           aes(label=gene_name), size=3, min.segment.length=0,
                           fontface="italic")+
  labs(x="PC3 gene loading", y="gene variance explained: slide",
       title="PRECAST (smoothed)")+
  theme_minimal()+theme(text=element_text(size=9))

p4<- ggplot(d1_se, aes(x=pc3, y=var_slide))+
  geom_point(size=.3, color="grey")+
  ggrepel::geom_text_repel(data=filter(d1_se, gene_name %in% genes2),
                           aes(label=gene_name), size=3, min.segment.length=0,
                           fontface="italic")+
  labs(x="PC3 gene loading", y="gene variance explained: slide",
       title="Seurat PC30")+
  theme_minimal()+theme(text=element_text(size=9))


# plots 5 and 6: pc3 variance explained and slide variance explained (all genes)
var1 = c("smoothed_k9_1663","sample_id","slide","seq")
d2_sm = tibble::rownames_to_column(as.data.frame(var.m_sm[,c("pc3",var1)]), var="gene_id") %>%
  tidyr::pivot_longer(cols=all_of(var1), names_to="expr_var", values_to="var_expl") %>%
  mutate(expr_var=factor(expr_var, levels=var1))
p6 <- ggplot(d2_sm, aes(x=pc3, y=var_expl))+
  geom_point(size=.1)+facet_wrap(vars(expr_var))+
  labs(x="PC3 variance explained", y="variance explained",
       title="PRECAST (smoothed) (all genes)")+
  theme_minimal()+theme(text=element_text(size=9))
var2 = c("seurat_label","sample_id","slide","seq")
d2_se = tibble::rownames_to_column(as.data.frame(var.m_se[,c("pc3",var2)]), var="gene_id") %>%
  tidyr::pivot_longer(cols=all_of(var2), names_to="expr_var", values_to="var_expl") %>%
  mutate(expr_var=factor(expr_var, levels=var2))
p7 <- ggplot(d2_se, aes(x=pc3, y=var_expl))+
  geom_point(size=.1)+facet_wrap(vars(expr_var))+
  labs(x="PC3 variance explained", y="variance explained",
       title="Seurat PC30 (all genes)")+
  theme_minimal()+theme(text=element_text(size=9))


# page 1 plots

g1 = arrangeGrob(p1, p2, p3, p4, rasterize(p6, dpi=150), rasterize(p7, dpi=150), ncol=2)




# slide vs pc3 pseudbulk sample embeddings
sort_slide = bind_rows(as.data.frame(colData(spe_sm))[,c("slide","pc3")],
          as.data.frame(colData(spe_se))[,c("slide","pc3")]) %>%
  group_by(slide) %>% summarise(avg_pc3=mean(pc3), med_pc3=median(pc3)) %>%
  arrange(med_pc3)

spe_sm$slide_f = factor(spe_sm$slide, levels=sort_slide$slide)
spe_se$slide_f = factor(spe_se$slide, levels=sort_slide$slide)

colpal = c("#b2df8a","#33a02c","#cab2d6","#6a3d9a")
names(colpal) = c("SC-TC","SKCCC","Psomagen-1","Psomagen-2")

p8.1 <- ggplot(as.data.frame(colData(spe_sm)), 
       aes(x=slide_f, y=pc3, color=smoothed_k9_1663))+
  geom_point()+scale_color_manual("PRECAST\n(smoothed)", values=cpList$smoothed.bright)+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

p8.2 <- ggplot(as.data.frame(colData(spe_se)), 
               aes(x=slide_f, y=pc3, color=seurat_label))+
  geom_point()+
  scale_color_manual("Seurat\nPC30", values=cpList$transfer.bright,
                     labels=c("M.V", names(cpList$transfer.bright)[2:8]))+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))


p9.1 <- ggplot(as.data.frame(colData(spe_sm)), 
                aes(x=slide_f, y=pc3, color=condition))+
  geom_point()+scale_color_manual("Dx", values=cpList$dx.pal)+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

p9.2 <- ggplot(as.data.frame(colData(spe_se)), 
                aes(x=slide_f, y=pc3, color=condition))+
  geom_point()+scale_color_manual("Dx", values=cpList$dx.pal)+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

# page 2 plots
g2 = arrangeGrob(p8.1, p8.2, p9.1, p9.2, ncol=2)


# seq and experimenter
#add vars of experimenter and date
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

#clean up experimenter
REDCap_MBv = REDCap_MBv[REDCap_MBv[["MBv_sample"]] %in% spe_sm$MBv_sample,]
REDCap_MBv$experimenter = gsub("Stephanie Page", "SCP", REDCap_MBv$experimenter)
REDCap_MBv$experimenter = gsub("Lina Oh", "Lina", REDCap_MBv$experimenter)
REDCap_MBv$experimenter = gsub("Svitlana Bach", "SB", REDCap_MBv$experimenter)

colData(spe_sm) = merge(colData(spe_sm), REDCap_MBv)
colData(spe_se) = merge(colData(spe_se), REDCap_MBv)


p10.1 <- ggplot(as.data.frame(colData(spe_sm)), 
               aes(x=slide_f, y=pc3, color=seq))+
  geom_point()+scale_color_manual("Seq. Core", values=colpal)+
  guides(colour = guide_legend(ncol = 2))+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

p10.2 <- ggplot(as.data.frame(colData(spe_se)), 
               aes(x=slide_f, y=pc3, color=seq))+
  geom_point()+scale_color_manual("Seq. Core", values=colpal)+
  guides(colour = guide_legend(ncol = 2))+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

p11.1 <- ggplot(as.data.frame(colData(spe_sm)), 
                aes(x=slide_f, y=pc3, color=experimenter))+
  geom_point()+scale_color_brewer(palette = "Dark2", labels=c("A","B","C","D"))+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

p11.2 <- ggplot(as.data.frame(colData(spe_se)), 
                aes(x=slide_f, y=pc3, color=experimenter))+
  geom_point()+scale_color_brewer(palette = "Dark2", labels=c("A","B","C","D"))+
  geom_boxplot(color="black", alpha=.5, outliers=F)+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=.5),
                   axis.title.x=element_blank(), legend.position="bottom",
                   legend.text=element_text(size=8), legend.title=element_text(size=9))

# page 3 plots
g3 = arrangeGrob(p10.1, p10.2, p11.1, p11.2, ncol=2)



# plots by date

col.date = c("#E2B1C6","#FF8A97","#E872AA","#941D61", #"oct 16","oct 23","nov 1","dec 18",
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
names(col.date) <- as.character(sort(unique(as.Date(spe_sm$date))))

p12.1 <- ggplot(as.data.frame(colData(spe_sm)), 
                aes(x=as.Date(date), y=pc3, color=date))+
  geom_point()+scale_color_manual(values=col.date, guide="none")+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank())
p12.2 <- ggplot(as.data.frame(colData(spe_se)), 
                aes(x=as.Date(date), y=pc3, color=date))+
  geom_point()+scale_color_manual(values=col.date, guide="none")+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank())

p13.1 <- ggplot(as.data.frame(colData(spe_sm)), 
                aes(x=as.Date(date), y=pc3, color=seq))+
  geom_point()+scale_color_manual("Seq. Core", values=colpal)+
  guides(colour = guide_legend(ncol = 2))+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank(), legend.position="bottom")
p13.2 <- ggplot(as.data.frame(colData(spe_se)), 
                aes(x=as.Date(date), y=pc3, color=seq))+
  geom_point()+scale_color_manual("Seq. Core", values=colpal)+
  guides(colour = guide_legend(ncol = 2))+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank(), legend.position="bottom")

p14.1 <- ggplot(as.data.frame(colData(spe_sm)), 
                aes(x=as.Date(date), y=pc3, color=experimenter))+
  geom_point()+scale_color_brewer(palette = "Dark2", labels=c("A","B","C","D"))+
  labs(y="pseudobulk sample PC3 embedding", title="PRECAST (smoothed)")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank(), legend.position="bottom")
p14.2 <- ggplot(as.data.frame(colData(spe_se)), 
                aes(x=as.Date(date), y=pc3, color=experimenter))+
  geom_point()+scale_color_brewer(palette = "Dark2", labels=c("A","B","C","D"))+
  labs(y="pseudobulk sample PC3 embedding", title="Seurat PC30")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),
                   axis.title.x=element_blank(), legend.position="bottom")

# page 4 plots
g4 = arrangeGrob(p12.1, p12.2, p13.1, p13.2, p14.1, p14.2, ncol=2)


# save all plots
pdf(file="plots/07_dx_DE/pc3_slide_covariate.pdf", height=11, width=8)
plot(g1)
plot(g2)
plot(g3)
plot(g4)
dev.off()

cat("\nPlots saved to: plots/07_dx_DE/pc3_slide_covariate.pdf\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
