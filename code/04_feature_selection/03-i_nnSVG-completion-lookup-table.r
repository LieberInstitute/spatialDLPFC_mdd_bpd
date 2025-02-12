logList = list.files("code/04_feature_selection/array_logs/")
lut = do.call(rbind, lapply(logList, function(x) {
  x1 = substr(x, start=27, stop=45)
  y = readLines(paste0("code/04_feature_selection/array_logs/",x), n=10)[10]
  z = readLines(paste0("code/04_feature_selection/array_logs/",x))
  y1 = substr(y, start=0, stop=13)
  error_codes=vector("list")
  if(length(grep("Error|error", z))>0) {
    if(length(grep("BiocParallel errors",z))>0) error_codes= c(error_codes, "BiocParallel")
    if(length(grep("BRISC",z))>0) error_codes= c(error_codes, "BRISC")
  }
  cbind.data.frame("log"=x1, "sample_id"=y1, 
                   "generate_weights"=file.exists(paste0("processed-data/04_feature_selection/per-sample_weights/",y1,"_spoon-weights.rda")),
                   "weighted_nnSVG"=file.exists(paste0("processed-data/04_feature_selection/per-sample_svgs/",y1,"_weighted-nnSVG-results.csv")),
                   "job_end"=length(grep("Job ends",z))>0,
                   "error"=paste(grep("Error|error", z),collapse="/"),
                   "error_code"=paste(error_codes, collapse="/"))
}))
write.table(lut, "processed-data/04_feature_selection/nnSVG-completion-lookup-table.txt", sep="\t", row.names = F)
