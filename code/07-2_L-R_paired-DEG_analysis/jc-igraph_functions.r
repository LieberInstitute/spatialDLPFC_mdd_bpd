getJC <- function(gmtSig_DFrame) {
  jcList = list()
  for(i in 1:nrow(gmtSig_DFrame)) {
    ref.set = gmtSig_DFrame[["leadingEdge"]][[i]]
    tmp.j = gmtSig_DFrame[-i,]
    out.mtx = matrix(NA, nrow=nrow(tmp.j), ncol=3, dimnames=list(tmp.j$pathway, c("size_query","size_overlap","jaccard")))
    for(j in 1:nrow(tmp.j)) {
      query.set = tmp.j[["leadingEdge"]][[j]]
      out.mtx[j, "size_query"]= length(query.set) 
      out.mtx[j, "size_overlap"] = length(intersect(query.set, ref.set))
      out.mtx[j, "jaccard"] = length(intersect(ref.set, query.set))/length(union(ref.set, query.set))
    }
    out.df = cbind.data.frame(as.data.frame(list("reference"=rep(gmtSig_DFrame[["pathway"]][i], nrow(out.mtx)),
                                                 "size_reference"=rep(length(ref.set), nrow(out.mtx)),
                                                 "query"=rownames(out.mtx))),
                              as.data.frame(out.mtx))
    rownames(out.df) <- NULL
    jcList[[gmtSig_DFrame[["pathway"]][i]]] = out.df
  }
  
  jc.df = do.call(rbind, jcList)
  rownames(jc.df) <- NULL
  return(jc.df)
}


pullJC <- function(gmtSet, .list_level) {
  #gmtSet reactList or WikiList
  #.list_level is either integer or one of c("padj","LR","paired")
  gmtPadj = lapply(gmtSet, pullDFrames, list_level=1)
  gmtPadj = setNames(do.call(c, gmtPadj), unlist(sapply(gmtPadj, names)))
  
  gmtLevel = lapply(gmtSet, pullDFrames, list_level=.list_level)
  gmtLevel = setNames(do.call(c, gmtLevel), unlist(sapply(gmtLevel, names)))
  
  #get jaccard coef for the leading edge gene sets
  do.call(rbind, lapply(names(gmtPadj), function(x) {
    getJC(gmtPadj[[x]]) %>% mutate(orig_group= x) %>%
      # filter to the sig gsea results that also contain a DEG at the specified list level
      filter(reference %in% gmtLevel[[x]]$pathway,
             query %in% gmtLevel[[x]]$pathway)
  })
  )
}

generateIGRAPH <- function(.orig_group, .gmt_db) {
  set.seed(123) #i think this is needed
  gmtSig_DFrame = gmtLR[[.orig_group]]
  filter_pathways = gmtPair[[.orig_group]]$pathway
  
  # node/ vertex data frame
  n.df = gmtSig_DFrame[,c("pathway","NES", grep("size_sig", colnames(gmtSig_DFrame), value=T))]
  colnames(n.df) = c("name","NES","size_sig")
  if(max(n.df$size_sig)<15) {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,max(n.df$size_sig)))
  } else {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,14))
  }
  n.df$term_level = factor(n.df$name %in% filter_pathways,
                           levels=c(FALSE, TRUE),
                           labels=c("LR","LR paired"))
  ## pretty names
  if(.gmt_db=="Reactome") {
    n.df$name2 = sapply(strsplit(n.df$name, " R-HSA-"), function(x) x[[1]])
  } 
  if(.gmt_db=="WikiPathways") {
    n.df$name2 = sapply(strsplit(n.df$name, " WP"), function(x) x[[1]])
  }
  n.df$name2 = sapply(n.df$name2, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  # edge data frame
  e.df= filter(jc.df, orig_group==.orig_group, jaccard>0) %>% select(from= reference, to= query, jc= jaccard)
  
  #make graph
  g <- graph_from_data_frame(e.df, directed=F, vertices=n.df)
  E(g)$weight <- e.df$jc
  
  V(g)$color <- as.character(factor(paste(sign(n.df$NES), n.df$term_level),
                                    levels=c("-1 LR","-1 LR paired","1 LR", "1 LR paired"),
                                    labels=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
  )
  V(g)$label.cex = .5
  V(g)$label.color = "black"
  V(g)$label = n.df$name2
  V(g)$size = n.df$size_scaled
  
  tmp = unlist(strsplit(.orig_group, split="_"))
  annot_tmp = as.character(factor(tmp[[1]], levels=c("sm","se"), labels=c("PRECAST (smoothed)","Seurat labels")))
  
  return(list(igraph= g,
              layout= layout_with_fr(g, weights=E(g)$jc*3),
              title_text= paste(annot_tmp, tmp[[2]]),
              sub_text= paste("Node size linearly corresponds to the number of", .orig_group,
                              "LR DEGs within leading edge.\nSmallest node =", 
                              min(n.df$size_sig), "DEG(s), largest node =", max(n.df$size_sig),"DEGs")
  ))
}
