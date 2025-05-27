library(SpatialExperiment)
library(dplyr)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")

avg.tbl = group_by(as.data.frame(colData(spe_pseudo)), precast_k9_1663, condition, sex) %>% summarise_at(c("nspots","sum","detected","subsets_mito_percent"), mean)

sd.tbl = group_by(as.data.frame(colData(spe_pseudo)), precast_k9_1663, condition, sex) %>% summarise_at(c("nspots","sum","detected","subsets_mito_percent"), sd)

n.tbl = group_by(as.data.frame(colData(spe_pseudo)), precast_k9_1663, condition, sex) %>% tally()

