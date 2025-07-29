setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(edgeR)
	library(BiocParallel)
})

set.seed(123)

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo)

spe_pseudo$pc3 = reducedDim(spe_pseudo)[,"PC3"]

cat("\nenrichment model: ~ seurat_label + condition + sex + nspots + pc3\n")
var_registration = "seurat_label"
covars = c("condition","sex","nspots","pc3")

#following guidance of spatialLIBD function to create contrast matrices
cluster_idx <- split(seq(along = spe_pseudo[[var_registration]]), spe_pseudo[[var_registration]])
modelList <- lapply(cluster_idx, function(x) {
            res <- rep(0, ncol(spe_pseudo))
            res[x] <- 1
            if (!is.null(covars)) {
                res_formula <-
                    eval(str2expression(paste(
                        "~", "res", "+", paste(covars, collapse = " + ")
                    )))
            } else {
                res_formula <- eval(str2expression(paste("~", "res")))
            }
            model.matrix(res_formula, data = colData(spe_pseudo))
})

#cat("\nvoomLmFit applied = TRUE\n")
#use bplapply to speed up
#fitList <- bplapply(modelList, function(x) {
#	voomLmFit(dge_pseudo, design=x, block=spe_pseudo$sample_id, save.plot=T)
#}, BPPARAM=MulticoreParam(workers=8))
#saveRDS(fitList, "processed-data/06_pseudobulk/voomLmFit-list_combat-seq_layer-enrichment_precast-k9-1663_covars-condition-sex-nspots.rda")
#cat("\nvoomLmFit objects with enrichment results saved to: processed-data/06_pseudobulk/voomLmFit-list_combat-seq_layer-enrichment_precast-k9-1663_covars-condition-sex-nspots.rda\n")


cat("\nvoom applied = FALSE\n")
cor_mod = model.matrix(~0 + seurat_label + condition + sex + nspots + pc3, data=colData(spe_pseudo))
corfit <- duplicateCorrelation(logcounts(spe_pseudo), design=cor_mod, block=spe_pseudo$sample_id)
fitList <- bplapply(modelList, function(x) {
	lmFit(logcounts(spe_pseudo), design=x, block=spe_pseudo$sample_id, correlation=corfit$consensus)
}, BPPARAM=MulticoreParam(workers=8))
saveRDS(fitList, "processed-data/06_pseudobulk/Seurat/lmFit-list_seurat-pc30-no-lowUMI_covars-condition-sex-nspots-pc3.rda")
cat("\nlmFit objects with enrichment results saved to: processed-data/06_pseudobulk/Seurat/lmFit-list_seurat-pc30-no-lowUMI_covars-condition-sex-nspots-pc3.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
