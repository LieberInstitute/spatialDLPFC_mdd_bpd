### modify enrichR .read_gmt to return it gmt file of enrichr database as a gson 
#https://github.com/wjawaid/enrichR/blob/master/R/functions.R#L232
.read_gmt <- function(db) {
  dbs <- listEnrichrDbs()
  if(!db %in% dbs$libraryName) stop("Requested database not in enrichR. Use enrichR::listEnrichrDbs() to check for available resources.")
  gmtDir = "code/enrichR_gmts"
  if(paste0(db,".rda") %in% list.files(gmtDir)) {
    gmt = readRDS(paste0(gmtDir, "/", db,".rda"))
  } else {
    base.address <- getOption("enrichR.base.address")
    url <- paste0(base.address, "geneSetLibrary?mode=text&libraryName=", db)
    tf <- tempfile(pattern = db, fileext = ".gmt")
    cat("   - Download GMT file...\n")
    tryCatch(download.file(url, tf, mode = "w", quiet = TRUE), 
             warning = function(warn) { message(warn); message("") },
             error = function(err) { message(err); message("") })
    gmt = read.gmt(tf)
    attr(gmt, 'snapshot') = c("retrieved"=format(Sys.time()), enrichR_version=packageVersion("enrichR"))
    attr(gmt, 'dbs_details') = as.list(dbs[grep(db, dbs$libraryName),])
    saveRDS(gmt, paste0(gmtDir, "/", db, ".rda"))
    cat("   - Saved to:", paste0(gmtDir, "/", db, ".rda"),"\n")
  }
  
  return(gmt)
}

### run fgsea for DE model results filtered to 1 dx*sex group and 1 cluster
runFGSEA <- function(filtered_df, .gmtList, adjp_threshold=.05, seed=123, ...) {
  set.seed(seed)
  cat("\nSeed check:", seed.check(),"\n")
  lf.df = group_by(filtered_df, gene_name) %>% summarise(avg_lf=mean(logFC)) %>% 
    arrange(desc(avg_lf))
  lf.order = lf.df$avg_lf
  names(lf.order) = lf.df$gene_name
  results = fgseaMultilevel(.gmtList, stats=lf.order, scoreType="std", minSize=10, maxSize=500, ...)
  set.seed(seed)
  sig = filter(results, padj<adjp_threshold)
  sig$leadingEdge2 = sapply(sig$leadingEdge, paste, collapse="/")
  cat("Seed check:", seed.check(),"\n")
  return(sig)
}
