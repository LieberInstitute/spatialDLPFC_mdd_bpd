library(SpatialExperiment)
library(HDF5Array)
library(raster)
library(escheR)
library(ggspavis)
library(ggplot2)
library(scuttle)
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
cdata2 = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")

identical(rownames(cdata), colnames(spe))

spe$genes_3MAD.outlier_binary = cdata$genes_3MAD.outlier_binary
spe$lg10.genes = log10(spe$sum_gene)
spe$umi100 = spe$sum_umi<=100
spe$true_edges = cdata2$true_edges
spe$problem_areas_binary = cdata2$problem_areas_binary
spe$plot_problem.areas = spe$true_edges | spe$problem_areas_binary


## poor quality section with problem area: 
#example_section = "V13B23-328_A1"
#example_section = "V13B23-302_B1"
#example_section = "V13Y10-021_A1"
#example_section = "V13B23-403_B1"

## good quality section with edge (and L1):
example_section = "V13B23-334_D1"



test = cdata[cdata$sample_id==example_section,c("array_row","array_col","genes_3MAD.outlier_binary")]
test2 = cdata[cdata$sample_id==example_section,c("array_row","array_col","in_tissue")]
spe_sub = spe[,spe$sample_id==example_section]
spe_sub = spe_sub[,spe_sub$in_tissue==T]

#normalization for MBP plotting
spe_sub <- computeLibraryFactors(spe_sub)
spe_sub <- logNormCounts(spe_sub)

plist = list()

#plot sum_genes
plist[[length(plist)+1]] = plotVisium(spe_sub, spots=TRUE, image=TRUE, facets=NULL, annotate="lg10.genes")+
    scale_fill_viridis_c("")+ggtitle("Detected genes (log10)")

#plot very low UMI
p <- plotVisium(spe_sub, spots=FALSE, image=TRUE, facets=NULL)+ggtitle("Sum UMI <= 100")
p1 = p %>% add_fill(var="umi100", point_size = 1.3)
plist[[length(plist)+1]] = p1+scale_fill_manual(values=c("grey","grey","red3"), guide="none")


#plot MBP expression
plist[[length(plist)+1]] = plotVisium(spe_sub, spots=TRUE, image=TRUE, facets=NULL, 
                        annotate=rownames(spe_sub)[rowData(spe_sub)$gene_name=="MBP"], assay = "logcounts", )+
  scale_fill_viridis_c("")+ggtitle("MBP expression (logcounts)")



#rasterize 3MAD outliers
colnames(test2) <- c("x","y","in_tissue")
odds = seq(1,max(test[,"array_col"]), by=2)

t1 = rasterFromXYZ(test)

cp1 = c("grey","black","white")
names(cp1) = c("0","1","off")
t1.df = raster::as.data.frame(t1, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
t1.df = merge(t1.df, test2)


#spot plot version
p <- plotVisium(spe_sub, spots=FALSE, image=TRUE, facets=NULL)+ggtitle("3MAD outlier: original")
p1 = p %>% add_fill(var="genes_3MAD.outlier_binary", point_size = 1.3)
plist[[length(plist)+1]] = p1+scale_fill_manual(values=c("grey","grey","black"), guide="none")


#focal transformation of 3MAD outliers
t2 = focal_transformations(t1)
t2.df = raster::as.data.frame(t2, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
t2.df = merge(t2.df, test2)


#spot plot version
t2.df2 = t2.df[,1:3]
colnames(t2.df2) <- c("array_row","array_col","transformed")
tmp = merge(colData(spe_sub), t2.df2)
rownames(tmp) = tmp$key
tmp$transformed = as.logical(tmp$transformed)
colData(spe_sub) = tmp[colnames(spe_sub),]

#have to remake visium because spe has new coldata
p <- plotVisium(spe_sub, spots=FALSE, image=TRUE, facets=NULL)+ggtitle("3MAD outlier: transformed")
p1 = p %>% add_fill(var="transformed", point_size = 1.3)
plist[[length(plist)+1]] = p1+scale_fill_manual(values=c("grey","grey","black"), guide="none")



#clump transformed 3MAD outliers
c1 = clump(t2, direction=8)
c1.df = raster::as.data.frame(c1, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
c1.df = merge(c1.df, test2)
#recode clump is so that they are more randomized
set.seed(123)
rnum = sample(1:max(c1.df$clumps, na.rm=T), max(c1.df$clumps, na.rm=T))
c1.df$clumps2 = factor(c1.df$clumps, levels=1:max(c1.df$clumps, na.rm=T), 
                       labels=paste(rnum, 1:max(c1.df$clumps, na.rm=T), sep="_"))

#spot plot version
c1.df2 = c1.df[,c("x","y","clumps2")]
colnames(c1.df2) <- c("array_row","array_col","clumps")
#necessary to jitter colors
c1.df2$clumps = ifelse(is.na(c1.df2$clumps), "0", as.character(c1.df2$clumps))
tmp = merge(colData(spe_sub), c1.df2)
rownames(tmp) = tmp$key
colData(spe_sub) = tmp[colnames(spe_sub),]

#have to remake visium because spe has new coldata
p <- plotVisium(spe_sub, spots=FALSE, image=TRUE, facets=NULL)+ggtitle("3MAD outlier: clumps")
p1 = p %>% add_fill(var="clumps", point_size = 1.3)
plist[[length(plist)+1]] = p1+scale_fill_manual(values=c("grey",viridisLite::turbo(n=max(c1.df$clumps, na.rm=T)),"grey"), guide="none")


#now highlight problem areas that were removed
p <- plotVisium(spe_sub, spots=FALSE, image=TRUE, facets=NULL)+ggtitle("3MAD outlier: problem areas")
p1 = p %>% add_fill(var="plot_problem.areas", point_size = 1.3)
plist[[length(plist)+1]] = p1+scale_fill_manual(values=c("grey","grey","black"), guide="none")

#plot lg10 genes without problem areas (same limits as before)
spe_sub$plot2 = ifelse(spe_sub$plot_problem.areas==T, NA, spe_sub$lg10.genes)
plist[[length(plist)+1]] <- plotVisium(spe_sub, spots=TRUE, image=TRUE, facets=NULL, annotate="plot2")+
  scale_fill_viridis_c("", limits=c(min(spe_sub$lg10.genes), max(spe_sub$lg10.genes)), na.value = "transparent")+
  ggtitle("Detected genes (log10) (problem areas removed)")



ggsave(file=paste0("plots/publication/example-",gsub("_", "-", example_section),"_3MAD-qc-approach_spot-plots.pdf"),
       gridExtra::marrangeGrob(plist, ncol=1, nrow=1, top=NULL),
       height=5, width=5)


### double check the precast results I ran when I excluded the flagged spots
library(Seurat)
load("processed-data/05_clustering/PRECAST/srt_no-problem-areas_precast_k-9_n1663.Rdata")
dim(seuInt@meta.data) #526409
sdata = seuInt@meta.data
rownames(sdata) = substr(rownames(sdata), start=15, stop=50)

prc.data = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
head(prc.data)
dim(prc.data) #535248
prc.data$no.problem.areas = NA
prc.data[rownames(sdata),"no.problem.areas"] = sdata$cluster
table(prc.data[,c("precast_k9_1663_f","no.problem.areas")])
#6 = GABA
#5 = L1
#7 = L2
#4 = L3.4 (1)
#9 = L3.4 (2)
#3 = L5
#1 = low UMI
#8 = WM

load("processed-data/05_clustering/PRECAST/srt_no-problem-areas-outliers_precast_k-9_n1663.Rdata")
dim(seuInt@meta.data) #523352
sdata = seuInt@meta.data
rownames(sdata) = substr(rownames(sdata), start=15, stop=50)

prc.data$no.problem.areas.outliers = NA
prc.data[rownames(sdata),"no.problem.areas.outliers"] = sdata$cluster
table(prc.data[,c("precast_k9_1663_f","no.problem.areas.outliers")])
#6 = L1
#2 = L2
#7 = L3.4
#1 = L5
#3 = L6
#5 = low UMI
#4 = Vasc
#8 = WM (1)
#9 = WM (2)

#### in both instances, there is still a low UMI cluster at k=9
