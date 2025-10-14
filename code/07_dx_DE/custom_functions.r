getTopTable <- function(eBayes_results, .coef, coef_name=NA) {
  if(.coef=="all") {
    t1 = topTable(eBayes_results, n=Inf)
  } else {
    t1 = topTable(eBayes_results, coef=.coef, n=Inf)
  }
  
  #t1$fdr = p.adjust(t1$P.Value, "fdr")
  if(!is.na(coef_name)) t1$coef = coef_name
  t1$gene_id = rownames(t1)
  t1$gene_name = rowData(spe_pseudo)[rownames(t1),"gene_name"]
  rownames(t1) <- NULL
  
  return(t1)
}


groupContrasts <- function(.fit_results, .comparisons, add_covar=NULL, return_contrasts=FALSE) {
  cont_mtx = matrix(0, nrow=ncol(coef(.fit_results)), ncol=length(.comparisons), 
                    dimnames = list(colnames(coef(.fit_results)), .comparisons))
  for(i in .comparisons) {
    seg = unlist(strsplit(i,"_"))
    if(length(seg)==1) {
      dx_refer = unlist(strsplit(seg, "\\."))[[1]]
      dx_compare = unlist(strsplit(seg, "\\."))[[2]]
      cont_mtx[grep(dx_refer, rownames(cont_mtx)),i] = -1
      cont_mtx[grep(dx_compare, rownames(cont_mtx)),i] = 1
    } else {
      dx = seg[grepl("\\.", seg)]
      dx_refer = unlist(strsplit(dx, "\\."))[[1]]
      dx_compare = unlist(strsplit(dx, "\\."))[[2]]
      clus = paste(seg[!grepl("\\.", seg)], collapse=".")
      #swap out dot placeholder
      clus = gsub("dot","\\.", clus)
      cont_mtx[grep(paste(dx_refer, clus, sep="."), rownames(cont_mtx)),i] = -1
      cont_mtx[grep(paste(dx_compare, clus, sep="."), rownames(cont_mtx)),i] = 1
    }
  }
  
  #add covar
  if(!is.null(add_covar)) {
    for(j in add_covar) {
      cont_mtx = cbind(cont_mtx, "tmp"=0)
      colnames(cont_mtx)[ncol(cont_mtx)] = j
      cont_mtx[j, j] = 1
    }
  }
  
  #change "dot" in rownames to match with coef names
  colnames(cont_mtx) <- gsub("dot","\\.",colnames(cont_mtx))
  
  if(return_contrasts) return(cont_mtx)
  contrasts.fit(.fit_results, cont_mtx)
  
}


sexTopTable <- function(.eBayes_results, phist=TRUE, cluster=NA, add_covar=NULL) {
  cnames = colnames(coef(.eBayes_results))
  if(!is.null(add_covar)) cnames= setdiff(cnames, add_covar)
  if(!is.na(cluster)) {cnames = grep(cluster, cnames, value=T)}
  resList = lapply(cnames, function(x) {
    t1 = getTopTable(.eBayes_results, x)
    t1$coef = x
    if(!is.na(cluster)) {
      t1$cluster = cluster
      t1$sex = factor(unlist(strsplit(x, "_"))[[2]], levels=c("F","M"))
      t1$group = factor(unlist(strsplit(x, "_"))[[3]], levels=c("NTC.MDD","NTC.BPD","MDD.BPD"))
    } else {
      t1$sex = factor(unlist(strsplit(x, "_"))[[1]], levels=c("F","M"))
      t1$group = factor(unlist(strsplit(x, "_"))[[2]], levels=c("NTC.MDD","NTC.BPD","MDD.BPD"))
    }
    return(t1)
  })
  
  tmp = do.call(rbind, resList)
  
  if(phist==FALSE) {
    return(tmp)
  } else {
    outList <- list("results"=tmp)
    outList$phist <- ggplot(tmp, aes(P.Value))+
      geom_histogram(bins=50)+facet_grid(cols=vars(sex), rows=vars(group))+
      theme_bw()
    return(outList)
  }
  
}
