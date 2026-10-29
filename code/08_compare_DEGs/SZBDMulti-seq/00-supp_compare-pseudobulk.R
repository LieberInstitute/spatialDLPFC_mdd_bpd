# this was run in interactive node
library(SingleCellExperiment)
library(dplyr)

# load csv export of raw counts for 3 sex chr genes and ncells for sanity check
df1 = read.csv("processed-data/08_compare_DEGs/check_sex-chrM_raw-counts.csv", row.names=1)

# load raw pseudobulk
load("processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res.Rdata")

genes1 = c("XIST","DDX3Y","RPS4Y1")
rowData(sce_pseudo)$gene_name = rowData(sce_pseudo)$featurekey
for(i in genes1) colData(sce_pseudo)[[i]] = counts(sce_pseudo)[rowData(sce_pseudo)$gene_name==i,]

cdata = as.data.frame(colData(sce_pseudo)) %>% group_by(individualID, Disorder, Biological_Sex) %>% summarise_at(c("ncells","XIST","DDX3Y","RPS4Y1"), sum)

# merge dfs, keep both ncells to check for identical
check.df = left_join(cdata, df1, by=c("individualID","Disorder","Biological_Sex"), suffix=c("_py","_seu"))

identical(check.df$ncells_py, check.df$ncells_seu) # TRUE
# so any mismatch is not the product of accidental nuclei filters

# look at plots
plot(log10(check.df$XIST+1), log10(check.df$XIST_counts+1), xlab="SCE counts", ylab="PY counts", main="XIST")
plot(log10(check.df$DDX3Y+1), log10(check.df$DDX3Y_counts+1), xlab="SCE counts", ylab="PY counts", main="DDX3Y")
plot(log10(check.df$RPS4Y1+1), log10(check.df$RPS4Y1_counts+1), xlab="SCE counts", ylab="PY counts", main="RPS4Y1")
# all are strongly correlated but not exactly identical so I need to go back through R pipeline and see what happened
