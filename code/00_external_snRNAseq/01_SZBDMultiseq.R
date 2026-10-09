### anndata::read_h5ad calls sub functions that create an specific type of virtual environment that can't be run in a batch session without additional tweaking
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(reticulate)
	library(anndata)
	library(Seurat)
})
# use sink to dump to log since in interactive mode
sink("code/00_external_snRNAseq/logs/SZBDMultiseq.log", append=FALSE, split=TRUE)

adata <- read_h5ad("raw-data/PsychEncode_backup/SZBDMulti-Seq_annotated.h5ad")
seu <- CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
# replace underscores with dashes (just like CreateSeuratObject did)
tmp = adata$var
rownames(tmp) = gsub("_","-",rownames(tmp))
seu[["RNA"]] <- AddMetaData(object = seu[["RNA"]], metadata = tmp)

# save space
rm(adata)

# collect all metadata
mdata = read.csv('https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt', sep ='\t')
mdata = mdata[mdata$Cohort=='SZBDMulti-Seq',]
colnames(mdata)[2] = "individualID"

s.mdata = seu@meta.data
s.mdata$key = rownames(s.mdata)

new.mdata = merge(s.mdata, mdata, sort=F)
rownames(new.mdata) <- new.mdata$key

new.mdata = new.mdata[rownames(s.mdata),]

# extract batch info from channel
new.mdata$Channel = unlist(lapply(strsplit(as.character(new.mdata$Channel), "_"), function(x) x[[1]]))
l1 = list(paste0("D19-", 7396:7403),
paste0("D19-", 6771:6778),
paste0("D19-", 5859:5866),
paste0("D19-", 7388:7395),
paste0("D19-", 7380:7387),
paste0("D19-", 6779:6786),
paste0("D19-", 4295:4302),
paste0("D19-", 7404:7411)
)
batch_encode = unlist(lapply(1:length(l1), function(x) rep(paste0("batch",x), length(l1[[x]]))))
channel_encode = unlist(l1)
new.mdata$batch = factor(new.mdata$Channel, levels=channel_encode, labels=batch_encode)
new.mdata = new.mdata[,c(1:7,25,8:24)]

# combine metadata
colSums(is.na(new.mdata))
stopifnot(identical(rownames(seu@meta.data), rownames(new.mdata)))
stopifnot(identical(seu@meta.data$nCount_RNA, new.mdata$nCount_RNA))

seu@meta.data <- new.mdata

# remove SZ
seu <- seu[,seu$Disorder %in% c("control","Bipolar Disorder")]
dim(seu) # 34291 312828


# add gene_name to row metadata
seu[["RNA"]]@meta.data$featurekey = rownames(seu)

# save object
saveRDS(seu, "processed-data/00_external_snRNAseq/seurat_SZBDMultiseq.rda")
write.csv(seu@meta.data, "processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_unfiltered_observations.csv", row.names=F)
write.csv(seu[["RNA"]]@meta.data, "processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_unfiltered_features.csv", row.names=F)

# nuclei filters
## numeric ages only
seu <- seu[,seu$Age_death!="89+"]
seu@meta.data$Age_death = as.numeric(seu@meta.data$Age_death)
dim(seu) # 34291 303771

## more than 1500 nuclei per sample
sort(table(seu$individualID))
### replicating original pipeline exactly, including mistaken inclusion of BD24 due to typo
seu <- seu[, !seu$individualID %in% c("CON14","CON24","BPD24","BD3")]
dim(seu) # 34291 301136

# filter to only genes with >100 nuclei according to initial metadata
seu = seu[seu[["RNA"]]@meta.data$n_cells>100,]
dim(seu) # 29887 301136

# save nuclei and gene filters
saveRDS(seu, "processed-data/00_external_snRNAseq/seurat_SZBDMultiseq_qc-filtered.rda")
write.csv(seu@meta.data, "processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_filtered_observations.csv", row.names=F)
write.csv(seu[["RNA"]]@meta.data, "processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_filtered_features.csv", row.names=F)

#split to process separately for cell type projection
seu_con = seu[,seu$Disorder=="control"]
dim(seu_con) # 29887 169536
save(seu_con, file="processed-data/00_external_snRNAseq/seurat_SZBDMultiseq_qc-filtered_control.Rdata")


seu_bd = seu[,seu$Disorder=="Bipolar Disorder"]
dim(seu_bd) # 29887 131600
save(seu_bd, file="processed-data/00_external_snRNAseq/seurat_SZBDMultiseq_qc-filtered_bipolar.Rdata")


cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessioninfo::session_info()
