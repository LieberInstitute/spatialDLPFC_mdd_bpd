library(Seurat)

## check for supp files with fetch_files = FALSE, then get with fetch_files = TRUE
#GEOquery::getGEOSuppFiles("GSE213982", fetch_files = FALSE)
## or download directly from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213982

# here we are loading files previously fetched that way
counts = Matrix::readMM("raw-data/Maitra-2023/GSE213982_combined_counts_matrix.mtx.gz")
obs = read.csv("raw-data/Maitra-2023/GSE213982_combined_counts_matrix_cells_columns.csv.gz")

# merge with metadata
df1 = as.data.frame(do.call(rbind, lapply(strsplit(obs$x, "\\."), function(x) x[1:3])))
colnames(df1) <- c("Sample", "barcode", "CellType")
df1$key = paste(gsub("-","_", df1$barcode), df1$Sample, sep="_")
rownames(df1) <- df1$key
colnames(counts) <- df1$key

mdata = read.csv("raw-data/Maitra-2023/Maitra-2023_ST1-sample-metadata.csv")
df1 = merge(df1, mdata, sort=F)
stopifnot(identical(colnames(counts), df1$key))
rownames(df1) <- df1$key


# feature data
var = read.csv("raw-data/Maitra-2023/GSE213982_combined_counts_matrix_genes_rows.csv.gz")

## use mbv rowData to match ENSEMBL to SYMBOL
mbv.var = read.csv("processed-data/93_globus/MBv_n119_features.csv.gz", row.names=1)
## remove chrM genes (already removed from Maitra 2023)
mbv.var = mbv.var[-grep("^MT-", mbv.var$gene_name),]
table(mbv.var$gene_name==var$x)
## now there are only 10 instances where genes names don't match and that is due to replicate ENSEMBL for multiple genes (Maitra gene names add '.1' to end)
var[var$x!=mbv.var$gene_name,]
mbv.var[mbv.var$gene_name!=var$x,"gene_name"]
## critically the order is identical so can use ensembl for rownames of counts
rownames(counts) <- mbv.var$gene_id

# make seurat object
seu <- CreateSeuratObject(counts = counts, meta.data = df1)
dim(seu) # 36588 160711

var = as.data.frame(mbv.var[,c("gene_id","gene_name")])
nuclei_more0 = rowSums(seu[["RNA"]]$counts>0)
var$n_nuclei = nuclei_more0

seu[["RNA"]] <- AddMetaData(object = seu[["RNA"]], metadata = var)

saveRDS(seu, "processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered.rda")

write.csv(seu@meta.data, "processed-data/00_external_snRNAseq/Maitra2023_control-MDD_filtered_observations.csv", row.names=F)
write.csv(seu[["RNA"]]@meta.data, "processed-data/00_external_snRNAseq/Maitra2023_control-MDD_unfiltered_features.csv", row.names=F)

# split by sex for label transfer since male and female samples were collected and processed years apart
## first filter to genes expressed in >50 nuclei (equivalent to SZBDMulti-seq 100 nuclei filter)
keep.genes = seu[["RNA"]]@meta.data$n_nuclei>50
seu = seu[keep.genes,]
dim(seu) # 28629 160711

write.csv(seu[["RNA"]]@meta.data, "processed-data/00_external_snRNAseq/Maitra2023_control-MDD_filtered_features.csv", row.names=F)

seu_m = seu[,seu$Sex=="Male"]
dim(seu_m) # 28629 79058
save(seu_m, file="processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered_male.Rdata")


seu_f = seu[,seu$Sex=="Female"]
dim(seu_f) # 28629 81653
save(seu_f, file="processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered_female.Rdata")
