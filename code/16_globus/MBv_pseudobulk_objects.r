setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
})

set.seed(123)

# domain-sp first
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
dim(spe_pseudo) # 21080   690

## reset rowData, add SVG and DE input columns
rowData(spe_pseudo) = rowData(spe_pseudo)[,1:7]

geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
rowData(spe_pseudo)$SVG = rowData(spe_pseudo)$gene_name %in% geneList$qual_genes

keep.sample <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
keep.domain <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
rowData(spe_pseudo)$DE_input = keep.sample==T & keep.domain==T

## update language
spe_pseudo$diagnosis = gsub("BPD", "BD", spe_pseudo$condition)
spe_pseudo$domain_sp = as.character(factor(spe_pseudo$smoothed_k9_1663, levels=c("L1","L2","L3.4","L5","L6","WM"), labels=c("L1","L2","L3/4","L5","L6","WM")))

colData(spe_pseudo) <- colData(spe_pseudo)[,c("sample_id", "brnum", "MBv_sample", "diagnosis", "sex", "age", "seq",
			"domain_sp", "nspots", "sum", "detected", "subsets_mito_percent")]

## round age to tenths for privacy
colData(spe_pseudo)$age = signif(colData(spe_pseudo)$age, 3)

saveRDS(spe_pseudo, "processed-data/16_globus/MBv_pseudobulk_domain-sp.rds")


# domain-ct next
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo) # 21175   942

## reset rowData, add SVG and DE input columns
rowData(spe_pseudo) = rowData(spe_pseudo)[,1:7]

rowData(spe_pseudo)$SVG = rowData(spe_pseudo)$gene_name %in% geneList$qual_genes

keep.sample <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
keep.domain <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
rowData(spe_pseudo)$DE_input = keep.sample==T & keep.domain==T

## update language
spe_pseudo$diagnosis = gsub("BPD", "BD", spe_pseudo$condition)
spe_pseudo$domain_ct = as.character(factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"), 
	labels=c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")))

colData(spe_pseudo) <- colData(spe_pseudo)[,c("sample_id", "brnum", "MBv_sample", "diagnosis", "sex", "age", "seq",
                        "domain_ct", "nspots", "sum", "detected", "subsets_mito_percent")]

## round age to tenths for privacy
colData(spe_pseudo)$age = signif(colData(spe_pseudo)$age, 3)

saveRDS(spe_pseudo, "processed-data/16_globus/MBv_pseudobulk_domain-ct.rds")


# domain-ct for eQTL next
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
dim(spe_pseudo) # 21174   943

## reset rowData, add SVG and DE input columns
rowData(spe_pseudo) = rowData(spe_pseudo)[,1:7]

rowData(spe_pseudo)$SVG = rowData(spe_pseudo)$gene_name %in% geneList$qual_genes

keep.sample <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
keep.domain <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
rowData(spe_pseudo)$GRNBOOST2_input = keep.sample==T & keep.domain==T

## update language
spe_pseudo$diagnosis = gsub("BPD", "BD", spe_pseudo$condition)
spe_pseudo$domain_ct = as.character(factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")))

colData(spe_pseudo) <- colData(spe_pseudo)[,c("sample_id", "brnum", "MBv_sample", "diagnosis", "sex", "age", "seq",
                        "domain_ct", "nspots", "sum", "detected", "subsets_mito_percent")]

## round age to tenths for privacy
colData(spe_pseudo)$age = signif(colData(spe_pseudo)$age, 3)

saveRDS(spe_pseudo, "processed-data/16_globus/MBv_pseudobulk_domain-ct-eQTL.rds")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
