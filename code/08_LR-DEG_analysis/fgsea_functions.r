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

### for use when distilling down to non-redundant terms
getTerms <- function(keyword, edges_DF) {
  keyword_terms = union(unique(grep(keyword, edges_DF$reference, ignore.case=T,value=T)),
                        unique(grep(keyword, edges_DF$query, ignore.case=T,value=T)))
  
  keyword_subset = filter(edges_DF, reference %in% keyword_terms | query %in% keyword_terms)
  #plot(density(keyword_subset[,"jaccard"]))
  keyword_subset2 = filter(keyword_subset, jaccard>.2)
  unique.terms = unique(c(keyword_subset2$reference, keyword_subset2$query, keyword_terms))
  unique.terms = unique(unlist(lapply(strsplit(unique.terms, " "), function(x) paste(x[-1], collapse=" "))))
  return(unique.terms)
}

### for use when combining annotations and removing redundant terms
mergeAnnotations <- function(.gmt_db, dx_sex_group) {
  source(paste0("code/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",dx_sex_group)), "_keylist.r"))
  #load smoothed 
  tmp2_sm = read.csv(paste0("processed-data/08_LR-DEG_analysis/", 
                            gsub("\\.","-", gsub("_","-",dx_sex_group)), "_smoothed-k9-1663_", .gmt_db, 
                            "_neuronal-fgsea_filtered-DEG-leadingeEdge.csv"))
  out2_sm = readRDS(paste0("processed-data/08_LR-DEG_analysis/", 
                           gsub("\\.","-", gsub("_","-",dx_sex_group)), "_smoothed-k9-1663_", .gmt_db, 
                           "_neuronal-fgsea_jc-igraph.rda"))
  keylist_sm = keylist[["smoothed"]][[.gmt_db]]
  keyterms_sm = lapply(keylist_sm$keywords, getTerms, edges_DF=out2_sm$edges)
  names(keyterms_sm) = keylist_sm$keywords
  names(keylist_sm$decidedterms) <- keylist_sm$keywords
  
  k_sm = do.call(rbind, lapply(keylist_sm$keywords, function(x) 
    data.frame("keyword"=rep(x, length(keyterms_sm[[x]])), 
               "term"=keyterms_sm[[x]], "keep"= keyterms_sm[[x]] %in% keylist_sm$decidedterms[[x]]))
  )
  #load seurat
  tmp2_se = read.csv(paste0("processed-data/08_LR-DEG_analysis/", 
                            gsub("\\.","-", gsub("_","-",dx_sex_group)), "_seurat-pc30_", .gmt_db, 
                            "_neuronal-fgsea_filtered-DEG-leadingeEdge.csv"))
  out2_se = readRDS(paste0("processed-data/08_LR-DEG_analysis/", 
                           gsub("\\.","-", gsub("_","-",dx_sex_group)), "_seurat-pc30_", .gmt_db, 
                           "_neuronal-fgsea_jc-igraph.rda"))
  keylist_se = keylist[["seurat"]][[.gmt_db]]
  keyterms_se = lapply(keylist_se$keywords, getTerms, edges_DF=out2_se$edges)
  names(keyterms_se) = keylist_se$keywords
  names(keylist_se$decidedterms) <- keylist_se$keywords
  
  k_se = do.call(rbind, lapply(keylist_se$keywords, function(x) 
    data.frame("keyword"=rep(x, length(keyterms_se[[x]])), 
               "term"=keyterms_se[[x]], "keep"= keyterms_se[[x]] %in% keylist_se$decidedterms[[x]]))
  )
  
  k_all = rbind(k_sm, k_se)
  
  #redundant terms
  redundant.terms = filter(k_all, keep==F)$term
  keep.terms = filter(k_all, keep==T)$term
  #make sure not accidentally removing anything i want to keep
  redundant.terms = setdiff(redundant.terms,keep.terms)
  
  #extra terms
  extra.terms = union(keylist_se$extra.terms, keylist_sm$extra.terms)
  #make sure not accidentally removing anything important
  extra.terms = setdiff(extra.terms, keep.terms)
  redundant.terms = c(redundant.terms, extra.terms)
  
  tmp2 = rbind(tmp2_sm, tmp2_se)
  tmp2$revised = !tmp2$term %in% redundant.terms
  return(tmp2)
}

