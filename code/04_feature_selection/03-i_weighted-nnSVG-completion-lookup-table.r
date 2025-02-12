logList = list.files("code/04_feature_selection/array_logs/")
lut = do.call(rbind, lapply(logList, function(x) {
  x1 = substr(x, start=27, stop=45)
  top = readLines(paste0("code/04_feature_selection/array_logs/",x), n=10)
  y = top[10]
  node = substr(top[7], start=10, stop=40)
  z = readLines(paste0("code/04_feature_selection/array_logs/",x))
  y1 = substr(y, start=0, stop=13)
  error_codes=vector("list")
  if(length(grep("Error|error", z))>0) {
    if(length(grep("BiocParallel errors",z))>0) error_codes= c(error_codes, "BiocParallel")
    if(length(grep("Error in BRISC",z))>0) error_codes= c(error_codes, "BRISC")
    if(length(grep("CANCELLED",z))>0) error_codes=c(error_codes, "cancelled")
    if(grep("Error|error",z)[1]==32) error_codes=c(error_codes,"JT error")
  }
  cbind.data.frame(#"spe_file"=y,
                   "log"=x1, "sample_id"=y1, "nodes"= node,
                   "generate_weights"=file.exists(paste0("processed-data/04_feature_selection/per-sample_weights/",y1,"_spoon-weights.rda")),
                   "weighted_nnSVG"=file.exists(paste0("processed-data/04_feature_selection/per-sample_svgs/",y1,"_weighted-nnSVG-results.csv")),
                   "job_end"=length(grep("Job ends",z))>0,
                   "error"=paste(grep("Error|error", z),collapse="/"),
                   "error_code"=paste(error_codes, collapse="/"))
}))

no_job_end = lut[lut$job_end==F,"sample_id"]
no_nnSVG_results = lut[lut$weighted_nnSVG==F,"sample_id"]
no_weights = lut[lut$generate_weights==F,"sample_id"]

#no_nnSVG_results will have everything but i just wanted to be safe
rerun_any = paste0(unique(c(no_job_end, no_nnSVG_results, no_weights)),".Rdata")

#write re-run file list
fileConn<-file("processed-data/04_feature_selection/per-sample_spe_RERUN_list.txt")
writeLines(rerun_any, fileConn)
close(fileConn)
