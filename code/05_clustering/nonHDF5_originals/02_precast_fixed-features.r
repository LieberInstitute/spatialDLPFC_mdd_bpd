setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
	library(here)
})
set.seed(123)

load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
##srt.sets = srt.sets[c("Br6529","Br6436","Br8221","Br8600","Br5799","Br6192")]
##srt.sets = srt.sets[c("Br5366","Br5916","Br6231","Br5639","Br5599",'Br5694')]
##feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))

###SVG based feature list
l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[
        unlist(lapply(l1, function(x) {
                if(dir.exists(here("processed-data","04_preprocessing",x))) return(FALSE)
                else {
                        split_1 = unlist(strsplit(x, split="_"))[[1]]
                        split.2 = unlist(strsplit(x, split="\\."))[[2]]
                        if(split_1=="nnSVG" & split.2=="csv") return(TRUE)
                        else {return(FALSE)}
                }
        }))
]
svg.df = do.call(rbind, lapply(l1, function(x) mutate(read.csv(here("processed-data","04_preprocessing",x)), file=x) %>% filter(padj<.05)))
svgs = unique(svg.df$gene_id)
cat("\nStarting SVGs (per-slide nnSVG analysis, p adj<.05 in any slide):",length(svgs),"genes\n")
#remove.svgs = distinct(svg.df, file, gene_id, gene_name) %>% group_by(gene_name, gene_id) %>% tally() %>% filter(n==1) %>% pull(gene_id)
#svg.more1 = setdiff(svg.nobias, remove.svgs)

###layer marker based feature list
#layer.markers.100 = read.csv(here("processed-data","04_preprocessing","EXT_TableS8_sig_genes_FDR5perc_enrichment.csv")) %>%
#	filter(spatial_domain_resolution=="Sp09") %>% group_by(test) %>% slice_min(n=100, fdr)
#lm100 = unique(layer.markers.100$ensembl)
#cat("\nStarting layer markers (Huuki-Meyers top 100 fdr per cluster/test):",length(lm100),"genes\n")

###exclude fixed sets of genes based on technical noise
bad.gene.list = readRDS(here("processed-data","04_preprocessing","biased_excluded_feature-list.rda"))
#svgs = setdiff(svgs, unlist(bad.gene.list$exclude[c("MTRN","DNAJ","HSP")]))
#cat("SVGs after excluding MTRN, DNAJ, and HSP genes:",length(svgs),"genes\n\n")
svgs = setdiff(svgs, unlist(bad.gene.list$exclude))
cat("SVGs after excluding MTRN, DNAJ & HSP, and RPS & RPL genes:",length(svgs),"genes\n\n")

###supplement SVGs with layer markers
#layer.markers = read.csv(here("processed-data","04_preprocessing","EXT_TableS8_sig_genes_FDR5perc_enrichment.csv")) %>%
#	filter(fdr<.0001, stat>0, spatial_domain_resolution=="Sp09") %>%
#	mutate(domain_simple=factor(test, 
#		levels=paste0("Sp09D0",c(1,2,3,5,8,4,7,6,9)), 
#		labels=c("L1","L1","L2","L3","L4","L5","L6","WM","WM")))
#lm = distinct(layer.markers, ensembl, domain_simple) %>% group_by(ensembl) %>% add_tally() %>%
#	filter(n==1, domain_simple!="L1", domain_simple!="WM") %>% pull(ensembl) %>% unique()
#cat("\nSupplementing with H-M layer markers (fdr<.0001, stat>0): sig. only in 1 domain, exclude L1 and WM\n")
#layer.markers2 = read.csv(here("processed-data","04_preprocessing","EXT_layer-markers_tableS5.csv"))
#lm2 = unique(filter(layer.markers2, fdr<.05)$ensembl)
#cat("\nSupplementing with Maynard S5 layer markers (fdr<.05)\n")
#cat("*** n supp layer markers total =",length(lm2),"\n")
#cat("*** n supp layer markers not already SVG =", length(setdiff(lm2, svgs)),"\n")
#svgs = union(svgs, lm)
#cat("Revised # SVGs with supplemented layer markers:",length(svgs),"\n")

biased.genes = readRDS(here("processed-data","04_preprocessing","bindev_biased-feature_5sd_list.rda"))
#cat("\nRemoving all 5SD biased genes:",length(unique(unlist(biased.genes))),"genes\n")
#svgs = setdiff(svgs, unique(unlist(biased.genes)))
#lm100 = setdiff(lm100, unique(unlist(biased.genes)))
#cat("\nRemoving only per-slide 5SD biased genes:",length(biased.genes[['nSD5.pooled.outlier']]),"genes\n")
#lm100 = setdiff(lm100, biased.genes[['nSD5.pooled.outlier']])
cat("\nRemoving only 5SD biased genes present in both lists:",length(intersect(biased.genes[[1]],biased.genes[[2]])),"genes\n")
svgs = setdiff(svgs, intersect(biased.genes[[1]],biased.genes[[2]]))
cat("Final SVG feature list after biased gene removal:",length(svgs),"genes\n\n")
#cat("Final layer marker feature list after biased gene removal:",length(lm100),"genes\n\n")

preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=svgs,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj <- PRECAST(PRECASTObj, K=7)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(PRECASTObj, file=here("processed-data","05_clustering","srt_precast-list_k-8-max-it-50_24samp_bindev.3k.Rdata"))
save(seuInt,file=here("processed-data","05_clustering","srt_precast_24samp_k-7_svgs-exclude-all-no-bias-both.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
