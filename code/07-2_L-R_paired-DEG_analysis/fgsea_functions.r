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


### perform fgsea on L-R results and filter based on L-R DEG lists
runFGSEA <- function(orig_group, match_group, gmtList, source=c("smoothed","seurat")) {
  set.seed(123) #before starting, make sure seed set
  #pass either sm_orig or se_orig to orig_group
  #gmtList is react.list or another gmt-based list
  if(source=="smoothed") {
    lf.df = filter(restr.results_sm, smoothed==orig_group, group==target_group, sex==target_sex) %>% 
      group_by(gene_name) %>% summarise(avg_lf=mean(logFC)) %>% arrange(desc(avg_lf))
    lf.order = lf.df$avg_lf
    names(lf.order) = lf.df$gene_name
  }
  if(source=="seurat") {
    lf.df = filter(restr.results_se, seurat_label_f==orig_group, group==target_group, sex==target_sex) %>%
      group_by(gene_name) %>% summarise(avg_lf=mean(logFC)) %>% arrange(desc(avg_lf))
    lf.order = lf.df$avg_lf
    names(lf.order) = lf.df$gene_name
  }
  
  group_dn = paste(target_group, target_sex, "dn", sep="_")
  group_up = paste(target_group, target_sex, "up", sep="_")
  
  sig.both.dir = c(pairList[[group_dn]][["LR_sig"]][[source]][[orig_group]],
                   pairList[[group_up]][["LR_sig"]][[source]][[orig_group]])
  paired.both.dir = c(pairList[[group_dn]][["LR_sig_paired"]][[match_group]], 
                      pairList[[group_up]][["LR_sig_paired"]][[match_group]])
  
  # fgsea
  results = fgseaMultilevel(gmtList, stats=lf.order, scoreType="std", minSize=10, maxSize=500)
  set.seed(123) #reset seed every time
  sig = filter(results, padj<.05)
  cat(nrow(sig),"results with padj<.05\n")
  sig.list = tibble::deframe(sig[,c("pathway","leadingEdge")])
  
  #filter to sig.any
  sig$leadingEdge2 = sapply(sig$leadingEdge, paste, collapse="/")
  sig_long = tidyr::separate_rows(sig, leadingEdge2, sep="/")
  suppressMessages({
    sig.filt_LR = mutate(sig_long, sig.LR= leadingEdge2 %in% sig.both.dir) %>%
      group_by(pathway, padj, log2err, NES, size, sig.LR) %>%
      summarise(leadingEdge2= paste(leadingEdge2, collapse="/"),
                size_sig.LR= sum(sig.LR)) %>%
      filter(size_sig.LR>0)
  })
  cat(nrow(sig.filt_LR),"results with padj<.05 and L-R DEG\n")
  #filter to sig.paired
  suppressMessages({
    sig.filt_paired = mutate(sig_long, sig.paired= leadingEdge2 %in% paired.both.dir) %>%
      group_by(pathway, padj, log2err, NES, size, sig.paired) %>%
      summarise(leadingEdge2= paste(leadingEdge2, collapse="/"),
                size_sig.paired= sum(sig.paired),
      ) %>%
      filter(size_sig.paired>0)
  })
  cat(nrow(sig.filt_paired),"results with padj<.05 and paired L-R DEG\n")
  outList = list(sig, sig.filt_LR, sig.filt_paired)
  names(outList) <- c(paste(source, orig_group, "padj", sep="_"),
                      paste(source, orig_group, "LR", sep="_"),
                      paste(source, orig_group, "paired", sep="_"))
  return(outList)
}


### format L-R results to run fgsea on both annotations and return list output
generateFGSEA <- function(.match_group, .gmtList) {
  #match group is one of the paired group names
  
  sm_orig = sm_match[grep(.match_group, names(sm_match))]
  se_orig = se_match[grep(.match_group, names(se_match))]
  
  smList = list()
  for(i in sm_orig) {
    cat("\n",.match_group, "> PRECAST (smoothed) >", i,"\n")
    smList[[i]] <- runFGSEA(i, .match_group, .gmtList, source="smoothed")
  }
  names(smList) = paste("sm", names(smList), sep="_")
  
  seList = list()
  for(i in se_orig) {
    cat("\n",.match_group, "> Seurat labels >", i, "\n")
    seList[[i]] <- runFGSEA(i, .match_group, .gmtList, source="seurat")
  }
  names(seList) = paste("se", names(seList), sep="_")
  
  return(list("smoothed"=smList, "seurat"=seList))
}


### extract fgsea results dframes of selected filter from all spatial domains/clusters
pullDFrames <- function(.match_groupList, list_level=c("padj","LR","paired")) {
  #list_level can be integer or c("")
  if(is.numeric(list_level)) N=list_level
  if(is.character(list_level)) {
    if(!list_level %in% c("padj","LR","paired")) {
      stop("'list_level' must be an integer or one of 'padj', 'LR', or 'paired'.")
    }
    N = as.numeric(as.character(factor(list_level, levels=c("padj","LR","paired"), labels=c(1,2,3))))
  }
  c(lapply(.match_groupList[["smoothed"]], function(x) ungroup(x[[N]])),
    lapply(.match_groupList[["seurat"]], function(x) ungroup(x[[N]]))
  )
}

