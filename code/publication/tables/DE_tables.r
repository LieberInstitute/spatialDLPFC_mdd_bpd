setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(writexl)
})

for(results_set in c("smoothed-k9-1663","seurat-pc30")) {
  if(results_set=="smoothed-k9-1663") set_name="domain-SP"
  if(results_set=="seurat-pc30") set_name="domain-CT"

  la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  la.degs = la.degs[,-grep("df\\.prior", colnames(la.degs))]

  lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  lr.degs = lr.degs[,-grep("df\\.prior", colnames(lr.degs))]

  lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
	mutate(domain="all", model="whole-tissue", annotation=set_name)

  lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
	rename(domain=cluster) %>%
	mutate(model="domain-restricted", annotation=set_name)

  lat = lat[,colnames(lrt)]

  deg.list = list("whole-tissue_F-test"= la.degs, "whole-tissue_t-test"=lat,
	"domain-restricted_F-test"= lr.degs, "domain-restricted_t-test"=lrt)

  write_xlsx(deg.list, path = paste0("processed-data/publication/supp_tables/", set_name, "_DE-results.xlsx"))

#  write.csv(la.degs, paste0("processed-data/publication/supp_tables/", set_name, "_whole-tissue_F-test-sig.csv"), row.names=F) 
#  write.csv(lr.degs, paste0("processed-data/publication/supp_tables/", set_name, "_domain-restricted_F-test-sig.csv"), row.names=F)
#  write.csv(lat, paste0("processed-data/publication/supp_tables/", set_name, "_whole-tissue_t-test.csv"), row.names=F)
#  write.csv(lrt, paste0("processed-data/publication/supp_tables/", set_name, "_domain-restricted_t-test.csv"), row.names=F)

 cat("\nSaved multi-sheet excel spreadsheet to:", paste0("processed-data/publication/supp_tables/", set_name, "_DE-results.xlsx"),"\n")
}


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
