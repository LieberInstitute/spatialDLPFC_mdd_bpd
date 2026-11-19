# add a check for image files in case spe was constructed from csvs and not spaceranger output

library(ggspavis)
library(gridExtra)

spe <- loadHDF5SummarizedExperiment(dir=spe_dir, prefix=spe_prefix)
spe$lg10.umi = log10(spe$sum_umi+1)

qc.areas = read.csv("Documents/ChatGPT/MBv_pipeline_revision/nextflow_outputs/test/02_QC-annotation/colData_edges-outlier-areas.csv", row.names=1)
stopifnot(identical(colnames(spe), rownames(qc.areas)))

spe$edge_outlier = qc.areas$edge_outlier
spe$outlier_area = qc.areas$outlier_area
spe$check_area = as.character(
    factor(paste(spe$edge_outlier, spe$outlier_area),
    levels = c("TRUE TRUE", "TRUE FALSE", "FALSE TRUE", "FALSE FALSE"),
    labels = c("edge artifact", "edge artifact", "low UMI area", "none"))
)
spe$check_area = ifelse(spe$sum_umi<100 & spe$check_area=="none", "<100 UMI", spe$check_area)
spe$check_area <- factor(spe$check_area, levels=c("edge artifact", "low UMI area", "<100 UMI", "none"))

sids = unique(colData(spe)[spe$check_area %in% c("edge artifact", "low UMI area"), "sample_id"])

slist <- lapply(sids, function(x) spe[,spe$sample_id==x & spe$in_tissue])
names(slist) <- sids
plist1 <- lapply(slist, function(x) suppressWarnings(plotVisium(x, annotate="lg10.umi", 
    image=T, spots=T,  point_size= .5, pal="viridis")))
plist2 <- lapply(slist, function(x) plotVisium(x, annotate="check_area", 
    image=T, spots=T, point_size= .5,
    pal = c("edge artifact"="navy", "low UMI area"="#00a000", "<100 UMI"="red", "none"="grey")))
plist.sort = do.call(c, lapply(sids, function(x) list(plist1[[x]], plist2[[x]])))
ggsave(file="plots/03_QC/check_edges-outlier-areas.pdf", 
    marrangeGrob(grobs=plist.sort, nrow=1, ncol=2, top= quote(sids[[g]])),
    height=3, width= 7)

ggsave(file="Downloads/check_edges-outlier-areas.pdf", 
    marrangeGrob(grobs=plist.sort, nrow=1, ncol=2, top= quote(sids[[g]])),
    height=3, width= 7)
