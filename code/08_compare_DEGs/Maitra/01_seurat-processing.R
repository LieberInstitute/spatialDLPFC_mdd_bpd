suppressPackageStartupMessages({
	library(Seurat)
})

# load MBv genes for filtering
rdata = read.csv("processed-data/06_pseudobulk/PRECAST/pseudobulk-sample-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)

### process males ###
cat("\nPROCESSING MALE SUBSET - MAITRA 2023...\n")
load("processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered_male.Rdata")
dim(seu_m) # 28629 79058

# first filter to MBv genes so I have fewer gene_names that replicate
keep.genes = intersect(seu_m[["RNA"]]@meta.data$gene_id, rownames(rdata))
cat("\nNumber of MBv-filtered Maitra genes:", length(keep.genes), "\n")
seu_m = seu_m[keep.genes,]

# need to regenerate object with gene names as rownames to match with SZBDMulti-seq
counts = seu_m[["RNA"]]@layers$counts
colnames(counts) = colnames(seu_m)
cdata = seu_m@meta.data[,4:34]
## HSPA14 duplicate gene name
fdata = seu_m[["RNA"]]@meta.data
fdata[grep("HSPA14", fdata$gene_name),"gene_name"] = paste(fdata$gene_name[grep("HSPA14", fdata$gene_name)], fdata$gene_id[grep("HSPA14", fdata$gene_name)], sep=".")
stopifnot(identical(length(fdata$gene_name), length(unique(fdata$gene_name))))
rownames(fdata) = fdata$gene_name
rownames(counts) = fdata$gene_name

seu_m = CreateSeuratObject(counts=counts, meta.data=cdata)
seu_m[["RNA"]] <- AddMetaData(object = seu_m[["RNA"]], metadata = fdata)
dim(seu_m)

# now process with seurat
set.seed(123)
cat("\nPre-processdataset...\n")
format(Sys.time())
seu_m <- SCTransform(seu_m)
seu_m <- RunPCA(seu_m)
seu_m <- FindNeighbors(seu_m, dims = 1:20, k.param=50)
seu_m <- RunUMAP(seu_m, dims=1:20)

# find clusters to determine if any single donor contributes to individual clusters
seu_m <- FindClusters(seu_m, resolution=0.3)

save(seu_m, file="processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT.Rdata")
cat("\nSaved male seurat to: processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT.Rdata\n")



# process females
cat("\nPROCESSING FEMALE SUBSET - MAITRA 2023...\n")
load("processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered_female.Rdata")
dim(seu_f) # 28629 81653

# first filter to MBv genes so I have fewer gene_names that replicate
keep.genes = intersect(seu_f[["RNA"]]@meta.data$gene_id, rownames(rdata))
cat("\nNumber of MBv-filtered Maitra genes:", length(keep.genes), "\n")
seu_f = seu_f[keep.genes,]

# need to regenerate object with gene names as rownames to match with SZBDMulti-seq
counts = seu_f[["RNA"]]@layers$counts
colnames(counts) = colnames(seu_f)
cdata =	seu_f@meta.data[,4:34]
## HSPA14 duplicate gene name
fdata = seu_f[["RNA"]]@meta.data
fdata[grep("HSPA14", fdata$gene_name),"gene_name"] = paste(fdata$gene_name[grep("HSPA14", fdata$gene_name)], fdata$gene_id[grep("HSPA14", fdata$gene_name)], sep=".")
stopifnot(identical(length(fdata$gene_name), length(unique(fdata$gene_name))))
rownames(fdata)	= fdata$gene_name
rownames(counts) = fdata$gene_name

seu_f = CreateSeuratObject(counts=counts, meta.data=cdata)
seu_f[["RNA"]] <- AddMetaData(object = seu_f[["RNA"]], metadata	= fdata)
dim(seu_f)

# now process with seurat
set.seed(123)
cat("\nPre-processdataset...\n")
format(Sys.time())
seu_f <- SCTransform(seu_f)
seu_f <- RunPCA(seu_f)
seu_f <- FindNeighbors(seu_f, dims = 1:20, k.param=50)
seu_f <- RunUMAP(seu_f, dims=1:20)

# find clusters to determine if any single donor contributes to individual clusters
seu_f <- FindClusters(seu_f, resolution=0.3)

save(seu_f, file="processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT.Rdata")
cat("\nSaved female seurat to: processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT.Rdata\n")


## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
