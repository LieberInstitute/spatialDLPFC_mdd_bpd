logList = list.files("code/04_feature_selection/array_logs/")
logList3 = logList[grep("min-100",logList)]
lut3 = do.call(rbind, lapply(logList3, function(x) {
  #run info
  x1 = substr(x, start=27, stop=49)
  z = readLines(paste0("code/04_feature_selection/array_logs/",x))
  y = z[10]
  y1 = substr(y, start=0, stop=13)
  #log info
  node = substr(z[7], start=10, stop=45)
  ##dim info
  keep_rows = na.exclude(as.numeric(strsplit(z[29], " ")[[1]]))
  if(length(keep_rows)!=1) keep_rows= keep_rows[2]
  min_genes = as.numeric(strsplit(z[34], " ")[[1]][8])
  min_spots = as.numeric(strsplit(z[38], " ")[[1]][8])
  ##error info
  error_codes=vector("list")
  if(length(grep("Error|error", z))>0) {
    if(length(grep("BiocParallel errors",z))>0) error_codes= c(error_codes, "BiocParallel")
    if(length(grep("Error in BRISC",z))>0) error_codes= c(error_codes, "BRISC")
    if(length(grep("CANCELLED",z))>0) error_codes=c(error_codes, "cancelled")
  }
  cbind.data.frame(#"spe_file"=y,
    "log"=x1, "sample_id"=y1, 
    "n_rows"=keep_rows, "min.genes_per.spot"=min_genes, "min.spots_per.gene"=min_spots,
    #"generate_weights"=file.exists(paste0("processed-data/04_feature_selection/per-sample_weights/",y1,"_spoon-weights.rda")),
    "standard_nnSVG"=file.exists(paste0("processed-data/04_feature_selection/per-sample_svgs/",y1,"_nnSVG-results.csv")),
    "job_end"=length(grep("Job ends",z))>0,"nodes"= node,
    "error"=paste(grep("Error|error", z),collapse="/"),
    "error_code"=paste(error_codes, collapse="/"))
}))

write.table(lut3, "processed-data/04_feature_selection/nnSVG_initial-completion_table.txt", sep="\t", row.names = F)

lut3[lut3$error!="",] 
## 2 biocparallel errors, 2 jobs timed out and were cancelled
rerun = paste0(lut3[lut3$error!="","sample_id"],".Rdata")
fileConn<-file("processed-data/04_feature_selection/per-sample_spe_RERUN_list.txt")
writeLines(rerun, fileConn)
close(fileConn)
