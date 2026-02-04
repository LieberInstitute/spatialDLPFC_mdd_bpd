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

source("code/08_dx-sex_DEG_analysis/fgsea_functions.r")

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
cat("\n\n")

(res_file = "smoothed-k9-1663")


# set variables for FGSEA
#(gmt_db = "Reactome")
##(gmt_db = "WikiPathways")
#(gmt_db = "GO-BP")
#(gmt_db = "GO-CC")
gmt_dbl = c("Reactome")

# automated from here
lr.results <- read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", res_file, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

#subset to WM
wm.results = filter(lr.results, cluster=="WM")

for(gmt_db in gmt_dbl) {

cat(paste0("\n",gmt_db))

if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
#if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")
if(gmt_db=="GO-CC") gmt = .read_gmt("GO_Cellular_Component_2025")
if(gmt_db=="GO-BP") gmt = .read_gmt("GO_Biological_Process_2025")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
  tibble::deframe()

resList <- list()
for(group_i in levels(wm.results$group)) {
  grp.filt = filter(wm.results, group==group_i)
  for(sex_i in levels(wm.results$sex)) {
    sex.filt = filter(grp.filt, sex==sex_i)
    resList[[paste(group_i, sex_i, sep="_")]] = sex.filt
  }
}

l1 = lapply(names(resList), function(x) {
	cat(paste0("\n",x))
	runFGSEA(resList[[x]], gmt.list)
})
names(l1) <- names(resList)

saveRDS(l1, file=paste0("processed-data/08_dx-sex_DEG_analysis/", res_file, "_", gmt_db, "_WM_LR-fgsea-list.rda"))
cat("\nSaved FGSEA results to:", paste0("processed-data/08_dx-sex_DEG_analysis/", res_file, "_", gmt_db, "_WM_LR-fgsea-list.rda"), "\n")
cat("\n\n\n")
}

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
