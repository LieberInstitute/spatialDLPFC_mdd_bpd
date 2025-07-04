setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SingleCellExperiment)
        library(edgeR)
	library(BiocParallel)
})

set.seed(123)

load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata")
dim(sce_pseudo)

sce_pseudo$sex = factor(sce_pseudo$Biological_Sex, levels=c("female","male"), labels=c("F","M"))

cat("\nenrichment model: ~ seurat_low.res + sex + ncells\n")
var_registration = "seurat_low.res"
covars = c("sex","ncells")

#following guidance of spatialLIBD function to create contrast matrices
cluster_idx <- split(seq(along = sce_pseudo[[var_registration]]), sce_pseudo[[var_registration]])
modelList <- lapply(cluster_idx, function(x) {
            res <- rep(0, ncol(sce_pseudo))
            res[x] <- 1
            if (!is.null(covars)) {
                res_formula <-
                    eval(str2expression(paste(
                        "~", "res", "+", paste(covars, collapse = " + ")
                    )))
            } else {
                res_formula <- eval(str2expression(paste("~", "res")))
            }
            model.matrix(res_formula, data = colData(sce_pseudo))
})

#cat("\nvoomLmFit applied = TRUE\n")
#use bplapply to speed up
#fitList <- bplapply(modelList, function(x) {
#	voomLmFit(dge_pseudo, design=x, block=spe_pseudo$sample_id, save.plot=T)
#}, BPPARAM=MulticoreParam(workers=8))
#saveRDS(fitList, "processed-data/06_pseudobulk/voomLmFit-list_combat-seq_layer-enrichment_precast-k9-1663_covars-condition-sex-nspots.rda")
#cat("\nvoomLmFit objects with enrichment results saved to: processed-data/06_pseudobulk/voomLmFit-list_combat-seq_layer-enrichment_precast-k9-1663_covars-condition-sex-nspots.rda\n")


cat("\nvoom applied = FALSE\n")
cor_mod = model.matrix(~0 + seurat_low.res + sex + ncells, data=colData(sce_pseudo))
corfit <- duplicateCorrelation(logcounts(sce_pseudo), design=cor_mod, block=sce_pseudo$individualID)
fitList <- bplapply(modelList, function(x) {
	lmFit(logcounts(sce_pseudo), design=x, block=sce_pseudo$individualID, correlation=corfit$consensus)
}, BPPARAM=MulticoreParam(workers=8))
saveRDS(fitList, "processed-data/06_pseudobulk/SZBDMulti-seq/lmFit-list_control-low-res_covars-sex-ncells.rda")
cat("\nlmFit objects with enrichment results saved to: processed-data/06_pseudobulk/SZBDMulti-seq/lmFit-list_control-low-res_covars-sex-ncells.rda\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
