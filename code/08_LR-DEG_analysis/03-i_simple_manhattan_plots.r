library(dplyr)
library(ggplot2)
library(enrichR)

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/08_LR-DEG_analysis/fgsea_functions.r")
source("code/08_LR-DEG_analysis/manhattan_functions.r")

gmt_db = "WikiPathways"
res_file = "smoothed-k9-1663"
clust_levels = c("L1","L2","L3.4","L5","L6","WM")
names(clust_levels) = clust_levels

#need a specific term (in future will loop over many terms)
term1 = "Glucocorticoid Receptor Pathway WP2880"

#need whole gmt list
if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
  tibble::deframe()

#need fgsea results
fgsea.results = readRDS("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"))
#test with 1 cluster
gs1 = wiki_sm[["NTC.BPD_F"]][["L5"]]

#need results: filtered to 1 dx*sex group and 1 cluster
restr.results <- read.csv(paste0("processed-data/07_dx_DE/layer-restricted-age_", res_file, "_rev-gene-input_compiled-results.csv"), row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         cluster=factor(cluster, levels=clust_levels),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))
#test with 1 cluster
df1 = filter(restr.results, group=="NTC.BPD", sex=="F", cluster=="L5")

# format dataframes for plotting
test = formatManhattan(df1, gs1, gmt.list, term1)

# plot manhattan plot
plotManhattan(test$manhattan)

# plot volcano with genes labelled by leadingEdge and in term or not 
plotVolcano(test$volcano)
