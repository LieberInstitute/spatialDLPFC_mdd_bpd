### anndata::read_h5ad calls sub functions that create an specific type of virtual environment that can't be run in a batch session without additional tweaking
### so I just ran the following code in an interactive session requiring 350GB of memory
### (killed when req 150. canceled job and used seff and saw it used 300GB)

#setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
#suppressPackageStartupMessages({
	library(Seurat)
	library(reticulate)
	library(anndata)
#})
set.seed(123)

#pulled from here: https://github.com/satijalab/seurat/issues/9072
adata <- read_h5ad("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/05_clustering/Seurat/adata_SZBDMulti-Seq_filtered.h5")

seu <- CreateSeuratObject(counts = t(adata$layers[["counts"]]), meta.data = adata$obs)

#replace underscores with dashes (just like CreateSeuratObject did)
tmp = adata$var
rownames(tmp) = gsub("_","-",rownames(tmp))
seu[["RNA"]] <- AddMetaData(object = seu[["RNA"]], metadata = tmp)

#save space
rm(adata)

#split to process separately
seu_con = seu[,seu$Disorder=="control"]
dim(seu_con) #29887 169536
save(seu_con, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control.Rdata")

seu_bd = seu[,seu$Disorder=="Bipolar Disorder"]
dim(seu_bd) #29887 131600
save(seu_bd, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar.Rdata")

