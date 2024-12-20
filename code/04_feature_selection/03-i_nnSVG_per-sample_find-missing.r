######### performed in interactive session to determine which slides/samples to re-run
library(HDF5Array)
library(DelayedArray)
library(SpatialExperiment)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spe$dummy_slide = as.character(spe$slide)
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

out.dir = "processed-data/04_feature_selection/per-sample_svgs"
svg.list = list.files("processed-data/04_feature_selection/per-sample_svgs")
length(svg.list)

missing.files = setdiff(paste0(spe$sample_id, "_nnSVG-results.csv"), svg.list)

#find logs matching to missing files to see what the errors are
log.list = grep("^nnSVG_per-sample", list.files("code/04_feature_selection/multi-file_logs"), value=T)
names(log.list) <- log.list

key = sapply(log.list, function(x) {
  ch = readLines(paste0("code/04_feature_selection/multi-file_logs/",x), n=19)[[19]]
  substr(ch, start=6, stop=15)
})

sapply(missing.files, function(x) {
  names(key)[grep(substr(x, start=0, stop=10), key)]
})
#26 x4, 27 x3, 12 x2, 17 x2, 24

#errors
# 26: Stop worker failed with the error: wrong args for environment subassignment (6060) (BiocParallel error)
# 27: Stop worker failed with the error: wrong args for environment subassignment (6058)
# 12: Stop worker failed with the error: wrong args for environment subassignment (6058)
# 17: Stop worker failed with the error: wrong args for environment subassignment (6061)
# 24: Stop worker failed with the error: wrong args for environment subassignment (6061)

#per-sample_re-run_list.txt
sort(unique(sapply(missing.files, function(x) strsplit(x, split="_")[[1]][[1]])))
