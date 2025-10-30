setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(fgsea)
	library(enrichR)
})

set.seed(123)
seed.check <- function() sample(1:10, size=4)
cat("\nInitial seed check for set.seed(123):\n")
seed.check()

source("code/08_LR-DEG_analysis/fgsea_functions.r")

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
cat("\n\n")


# set variables for FGSEA
#(gmt_db = "Reactome")
#(gmt_db = "WikiPathways")
#(gmt_db = "GO-BP")
(gmt_db = "GO-CC")

(res_file = "smoothed-k9-1663")
(clust_levels = c("L1","L2","L3.4","L5","L6","WM"))

#(res_file = "seurat-pc30")
#(clust_levels = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))

names(clust_levels) = clust_levels

# automated from here
restr.results <- read.csv(paste0("processed-data/07_dx_DE/layer-restricted-age_", res_file, "_rev-gene-input_compiled-results.csv"), row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         cluster=factor(cluster, levels=clust_levels),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))


if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")
if(gmt_db=="GO-CC") gmt = .read_gmt("GO_Cellular_Component_2025")
if(gmt_db=="GO-BP") gmt = .read_gmt("GO_Biological_Process_2025")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
  tibble::deframe()

resList <- list()
for(group_i in levels(restr.results$group)) {
  grp.filt = filter(restr.results, group==group_i)
  for(sex_i in levels(restr.results$sex)) {
    sex.filt = filter(grp.filt, sex==sex_i)
    resList[[paste(group_i, sex_i, sep="_")]] = sex.filt
  }
}

l2 = lapply(names(resList), function(x) {
  cat(paste0("\n",x))
  l1 = lapply(clust_levels, function(y) {
    tmp = filter(resList[[x]], cluster==y)
    cat(paste0("\n",y))
    runFGSEA(tmp, gmt.list)
  })
  cat("\n\n\n")
  return(l1)
})
names(l2) <- names(resList)

saveRDS(l2, file=paste0("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"))
cat("\nSaved FGSEA results to:", paste0("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"), "\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
