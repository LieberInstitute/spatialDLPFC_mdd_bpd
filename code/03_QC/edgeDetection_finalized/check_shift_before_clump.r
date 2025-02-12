library(raster)
library(dplyr)
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)

#look at difference in edges with shift before clump
source("code/03_QC/edgeDetection_finalized/raster_edge_functions2.r")
sampleList = unique(cdata$sample_id)
names(sampleList) <- sampleList

genes_edges = lapply(sampleList, function(x) {
  tmp = cdata[cdata$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
  clumpEdges(tmp[,-1], rownames(tmp)[tmp$in_tissue==FALSE])
})
table(sapply(genes_edges, length)>0) 
sort(names(genes_edges)[sapply(genes_edges, length)>0])


genes_edges2 = lapply(sampleList, function(x) {
  tmp = cdata[cdata$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
  clumpEdges_shifted(tmp[,-1], rownames(tmp)[tmp$in_tissue==FALSE])
})
table(sapply(genes_edges2, length)>0) 
sort(names(genes_edges2)[sapply(genes_edges2, length)>0])

setdiff(names(genes_edges)[sapply(genes_edges, length)>0], names(genes_edges2)[sapply(genes_edges2, length)>0])
#"V13B23-302_C1"


genes_probs = lapply(sampleList, function(x) {
  tmp = cdata[cdata$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
  problemAreas(tmp[,-1], 
               rownames(tmp)[tmp$in_tissue==FALSE], uniqueIdentifier=x)
})
genes_probs = do.call(rbind, genes_probs)
cdata$problem_areas_genes.id2 = NA
cdata[genes_probs$spotcode,"problem_areas_genes.id2"] = genes_probs$clumpID
cdata$problem_areas_genes.size2 = 0
cdata[genes_probs$spotcode,"problem_areas_genes.size2"] = genes_probs$clumpSize



genes_probs2 = lapply(sampleList, function(x) {
  tmp = cdata[cdata$sample_id==x,c("in_tissue","array_row","array_col","genes_3MAD.outlier_binary")]
  problemAreas_shifted(tmp[,-1], 
               rownames(tmp)[tmp$in_tissue==FALSE], uniqueIdentifier=x)
})
genes_probs2 = do.call(rbind, genes_probs2)
cdata$problem_areas_genes.id3 = NA
cdata[genes_probs2$spotcode,"problem_areas_genes.id3"] = genes_probs2$clumpID
cdata$problem_areas_genes.size3 = 0
cdata[genes_probs2$spotcode,"problem_areas_genes.size3"] = genes_probs2$clumpSize



#doesn't change which spots were identified
group_by(cdata, sample_id) %>% summarise(spots_in_areas1 = sum(!is.na(problem_areas_genes.id2)),
                                         spots_in_areas2 = sum(!is.na(problem_areas_genes.id3))) %>% 
  filter(spots_in_areas1 != spots_in_areas2)

#just changes how they were clumped
filter(cdata, sample_id=="V13Y10-023_A1") %>% distinct(problem_areas_genes.id2, problem_areas_genes.size2) %>% slice_max(problem_areas_genes.size2, n=5)
filter(cdata, sample_id=="V13Y10-023_A1") %>% distinct(problem_areas_genes.id3, problem_areas_genes.size3) %>% slice_max(problem_areas_genes.size3, n=5)

#what about their prop
cdata$lowumi = cdata$sum_umi<=100

filter(cdata, in_tissue==T, sample_id=="V13Y10-023_A1") %>% group_by(problem_areas_genes.id2, problem_areas_genes.size2) %>% 
  summarise(n_spots_in_areas1=n(), prop_lowumi1=sum(lowumi)/n_spots_in_areas1) %>% ungroup() %>% slice_max(n_spots_in_areas1, n=5)
filter(cdata, in_tissue==T, sample_id=="V13Y10-023_A1") %>% group_by(problem_areas_genes.id3, problem_areas_genes.size3) %>% 
  summarise(n_spots_in_areas2=n(), prop_lowumi2=sum(lowumi)/n_spots_in_areas2) %>% ungroup() %>% slice_max(n_spots_in_areas2, n=5)
#shift before clump does expand size but in this case it results in area i want to be removed to now be <50% lowumi 
