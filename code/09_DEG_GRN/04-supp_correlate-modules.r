setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
  library(gridExtra)
})
set.seed(123)

# load in GRN adjacency output for correlations
lg.mask = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")

# load in module sets
#modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_interaction-modules.csv")
modules2 = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules.csv")

# empty input matrix
each_module = unique(modules2$TF)

input.mtx = matrix(NA, nrow=length(each_module), ncol=length(each_module), dimnames = list(each_module, each_module))

for(i in 1:(nrow(input.mtx)-1)) {
  complete.corr = colnames(input.mtx)[(i+1):ncol(input.mtx)]
  i=rownames(input.mtx)[i]
  for(j in complete.corr) {
    tmp1 = filter(modules2, TF==i)
    tmp2 = filter(modules2, TF==j)
    geneList = setdiff(union(tmp1$target, tmp2$target),
                       c(i,j))
    #pull the raw adj (before any filtering) and format for correlation
    check = filter(lg.mask, TF %in% c(i,j), target %in% geneList) %>%
      mutate(modules_key=factor(TF, levels=c(i,j), labels=c("I","J"))) %>%
      select(modules_key, target, importance) %>% tidyr::pivot_wider(names_from="modules_key", values_from="importance", values_fill=0)
    corrij= cor(check$I, check$J)
    input.mtx[i,j] = corrij
    input.mtx[j,i] = corrij
  }
}
write.csv(input.mtx, "processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv")
cat("\nSaved pairwise target correlation matrix for all DEG modules to: processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_DEG-modules_target-pearson-correlation.csv\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
