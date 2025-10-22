jcoef <- function(x, y) {
  length(intersect(x, y))/length(union(x, y))
}

jaccardFromList <- function(namedList) {
  d1 = data.frame("reference"=NA, "query"=NA, "size_reference"=NA, "size_query"=NA, "jaccard"=NA)
  nList = names(namedList)
  while(length(nList)>1) {
    ref = nList[1]
    query = nList[-1]
    d2 = do.call(rbind, lapply(query, function(x) {
      data.frame("reference"=ref, "query"=x, 
                 "size_reference"=length(namedList[[ref]]), "size_query"=length(namedList[[x]]),
                 "jaccard"=jcoef(namedList[[ref]], namedList[[x]]))
    }))
    d1 = rbind(d1, d2)
    nList = query
  }
  return(d1[-1,])
}

formatJaccardIGRAPH <- function(fgsea_results, .group, .cluster) {
  if(length(.cluster)>1) {
    return(formatJaccardIGRAPH_all(fgsea_results, .group, .cluster))
  }
  specific_fgsea_results = fgsea_results[[.group]][[.cluster]]
  #specific_fgsea_results is the list level of the fgsea results list for 1 dx*group and 1 cluster
  ledge = tibble::deframe(specific_fgsea_results[,c("pathway","leadingEdge")])
  jc.df = jaccardFromList(ledge)
  e.df= filter(jc.df, jaccard>0)

  
  lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]],
                  lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]])
  la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                  saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
  
  n.df = data.frame("term"=names(ledge),
                    "NES"=specific_fgsea_results$NES,
                    "LR_in_leadingEdge"=sapply(lapply(ledge, intersect, y=lr.both.dir), length), 
                    "LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=la.both.dir), length),
                    "LR.or.LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=union(la.both.dir, lr.both.dir)), length), 
                    "source"=rep(unlist(strsplit(res_file, "-"))[[1]], length(ledge)),
                    "cluster"=rep(y, length(ledge)),
                    row.names=NULL)
  ## pretty names
  if(gmt_db=="Reactome") {
    n.df$name = sapply(strsplit(n.df$term, " R-HSA-"), function(x) x[[1]])
  } 
  if(gmt_db=="WikiPathways") {
    n.df$name = sapply(strsplit(n.df$term, " WP"), function(x) x[[1]])
  }
  n.df$name = sapply(n.df$name, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  return(list("nodes"=n.df, "edges"=e.df))
}

### alt version for if multiple clusters are asked at once, it will calculate jaccard coeff of terms in different clusters
formatJaccardIGRAPH_all <- function(fgsea_results, .group, .clusters) {
  #use this version when .clusters is a character vector of more than 1 level and it will find the jaccard coef across terms of different levels
  specific_fgsea_results = lapply(fgsea_results[[.group]][.clusters], function(x) tibble::deframe(x[,c("pathway","leadingEdge")]))
  n2 = names(specific_fgsea_results)
  ledge = do.call(c, lapply(n2, function(x) {
    l1 = specific_fgsea_results[[x]]
    names(l1) <- paste(x, names(l1))
    return(l1)
    }))
  #specific_fgsea_results is the list level of the fgsea results list for 1 dx*group and 1 cluster
  jc.df = jaccardFromList(ledge)
  e.df= filter(jc.df, jaccard>0)
  
  merge.n.df = do.call(rbind, lapply(.clusters, function(y) {
    lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]],
                    lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[y]])
    la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                    saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
    
    n.df = data.frame("term"=paste(y, names(specific_fgsea_results[[y]])),
                      "NES"=fgsea_results[[.group]][[y]]$NES, 
                      "LR_in_leadingEdge"=sapply(lapply(specific_fgsea_results[[y]], intersect, y=lr.both.dir), length), 
                      "LA_in_leadingEdge"=sapply(lapply(specific_fgsea_results[[y]], intersect, y=la.both.dir), length),
                      "LR.or.LA_in_leadingEdge"=sapply(lapply(specific_fgsea_results[[y]], intersect, 
                                                              y=union(la.both.dir, lr.both.dir)), length), 
                      "source"=rep(unlist(strsplit(res_file, "-"))[[1]], length(specific_fgsea_results[[y]])),
                      "cluster"=rep(y, length(specific_fgsea_results[[y]])),
                      row.names=NULL)
    ## pretty names
    if(gmt_db=="Reactome") {
      n.df$name = sapply(strsplit(n.df$term, " R-HSA-"), function(x) x[[1]])
    } 
    if(gmt_db=="WikiPathways") {
      n.df$name = sapply(strsplit(n.df$term, " WP"), function(x) x[[1]])
    }
    n.df$name = sapply(n.df$name, function(x) paste(strwrap(x, width=20), collapse="\n"))
    return(n.df)
  }))
  
  return(list("nodes"=merge.n.df, "edges"=e.df))
}


### generate the igraph objects for plotting
generateIGRAPH <- function(formatted_jaccard) {
  if(length(unique(formatted_jaccard$nodes$cluster))>1) {
    return(generateIGRAPH_all(formatted_jaccard))
  }
  set.seed(123) # this is needed for consistent layouts
  # edge data frame
  e.df= select(formatted_jaccard$edges, from= reference, to= query, jc= jaccard)
  n.df = select(formatted_jaccard$nodes, name=term, label=name, NES, size_sig=LR.or.LA_in_leadingEdge)
  if(max(n.df$size_sig)<15) {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,max(n.df$size_sig)))
  } else {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,14))
  }
  
  #make graph
  g <- graph_from_data_frame(e.df, directed=F, vertices=n.df)
  E(g)$weight <- e.df$jc
  
  V(g)$color <- as.character(factor(paste(sign(n.df$NES), n.df$size_sig>0),
                                    levels=c("-1 FALSE", "-1 TRUE", "1 FALSE", "1 TRUE"),
                                    labels=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
  )
  V(g)$label.cex = .5
  V(g)$label.color = "black"
  V(g)$label = n.df$label
  V(g)$size = n.df$size_scaled
  annot_tmp = as.character(factor(unique(formatted_jaccard$nodes$source), 
                                  levels=c("smoothed","seurat"), 
                                  labels=c("PRECAST (smoothed)","Seurat labels")))
  
  outList = list(igraph= g,
                 layout= layout_with_fr(g, weights=E(g)$jc*3),
                 title_text= paste(annot_tmp, unique(formatted_jaccard$nodes$cluster)),
                 sub_text= paste("Node size linearly corresponds to the number of LR or LA DEGs within leading edge.\nSmallest node =", 
                                 min(n.df$size_sig), "DEG(s), largest node =", max(n.df$size_sig),"DEGs"))
  return(outList)
}

### alternative version for if terms of different clusters were checked against one another

generateIGRAPH_all <- function(formatted_jaccard) {
  set.seed(123) # this is needed for consistent layouts
  # edge data frame
  e.df= select(formatted_jaccard$edges, from= reference, to= query, jc= jaccard)
  n.df = select(formatted_jaccard$nodes, name=term, label=name, outline= cluster, NES, size_sig=LR.or.LA_in_leadingEdge)
  if(max(n.df$size_sig)<15) {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,max(n.df$size_sig)))
  } else {
    n.df$size_scaled= scales::rescale(n.df$size_sig, to=c(3,14))
  }
  
  #make graph
  g <- graph_from_data_frame(e.df, directed=F, vertices=n.df)
  E(g)$weight <- e.df$jc
  
  V(g)$color <- as.character(factor(paste(sign(n.df$NES), n.df$size_sig>0),
                                    levels=c("-1 FALSE", "-1 TRUE", "1 FALSE", "1 TRUE"),
                                    labels=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
  )
  V(g)$label.cex = .5
  V(g)$label.color = "black"
  V(g)$label = n.df$label
  V(g)$size = n.df$size_scaled
  V(g)$frame.color= as.character(factor(n.df$outline, levels=names(cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]]),
                                               labels=cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]]))
  annot_tmp = as.character(factor(unique(formatted_jaccard$nodes$source), 
                                  levels=c("smoothed","seurat"), 
                                  labels=c("PRECAST (smoothed)","Seurat labels")))
  
  outList = list(igraph= g,
                 layout= layout_with_fr(g, weights=E(g)$jc*3),
                 title_text= annot_tmp,
                 sub_text= paste("Node size linearly corresponds to the number of LR or LA DEGs within leading edge.\nSmallest node =", 
                                 min(n.df$size_sig), "DEG(s), largest node =", max(n.df$size_sig),"DEGs"))
  return(outList)
}
