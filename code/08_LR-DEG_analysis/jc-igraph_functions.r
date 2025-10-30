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

  
  lr.both.dir = c(lrList[[paste0(x,"_dn")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[.cluster]],
                  lrList[[paste0(x,"_up")]][["LR_sig"]][[unlist(strsplit(res_file, "-"))[[1]]]][[.cluster]])
  la.both.dir = c(saveList[["sig_genes"]][[paste0(x,"_dn")]][[paste0("adj_", substr(res_file, start=0, stop=2))]],
                  saveList[["sig_genes"]][[paste0(x,"_up")]][[paste0("adj_", substr(res_file, start=0, stop=2))]])
  
  n.df = data.frame("term"=names(ledge),
                    "NES"=specific_fgsea_results$NES,
                    "LR_in_leadingEdge"=sapply(lapply(ledge, intersect, y=lr.both.dir), length), 
                    "LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=la.both.dir), length),
                    "LR.or.LA_in_leadingEdge"=sapply(lapply(ledge, intersect, y=union(la.both.dir, lr.both.dir)), length), 
                    "source"=rep(unlist(strsplit(res_file, "-"))[[1]], length(ledge)),
                    "cluster"=rep(.cluster, length(ledge)),
                    row.names=NULL)
  ## pretty names
  if(gmt_db=="Reactome") {
    n.df$name = sapply(strsplit(n.df$term, " R-HSA-"), function(x) x[[1]])
  } 
  if(gmt_db=="WikiPathways") {
    n.df$name = sapply(strsplit(n.df$term, " WP"), function(x) x[[1]])
  }
  if(gmt_db %in% c("GO-BP","GO-CC")) {
    n.df$name = sapply(strsplit(n.df$term, " \\(GO:"), function(x) x[[1]])
  }
  n.df$name = sapply(n.df$name, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  return(list("nodes"=n.df, "edges"=e.df))
}

### alt version for if multiple clusters are asked at once, it will calculate jaccard coeff of terms in different clusters
formatJaccardIGRAPH_all <- function(fgsea_results, .group, .clusters) {
  #use this version when .clusters is a character vector of more than 1 level and it will find the jaccard coef across terms of different levels
  specific_fgsea_results = lapply(fgsea_results[[.group]][.clusters], function(x) tibble::deframe(x[,c("pathway","leadingEdge")]))
  #remove any empty list levels
  llength = sapply(specific_fgsea_results, length)
  .clusters = names(llength)[llength>0]
  specific_fgsea_results = specific_fgsea_results[.clusters]

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
    if(gmt_db %in% c("GO-BP","GO-CC")) {
      n.df$name = sapply(strsplit(n.df$term, " \\(GO:"), function(x) x[[1]])
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
  set.seed(123)
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
  set.seed(123)
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
  set.seed(123)
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
  set.seed(123)
  return(outList)
}


### different graph, with genes and terms
term2GeneIGRAPH <- function(unique.terms, leadingEdge_DF, res_file=res_file) {
  set.seed(123) 
  filtered_df = filter(leadingEdge_DF, term %in% unique.terms)
  
  #check that all in same direction
  st = table(sign(filtered_df$NES))
  if(dim(st)>1) warning("Not all terms in same direction (mix of depleted and enriched)")
  
  p.terms = tibble::deframe(distinct(filtered_df, term, name))[unique.terms]
  p.terms = sapply(p.terms, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  s.nes = tibble::deframe(distinct(filtered_df, term, sign(NES)))[unique.terms]
  
  n.df = data.frame("name"=c(unique.terms, unique(filtered_df$gene_name)),
                    "label"=c(p.terms, unique(filtered_df$gene_name)),
                    "type"=c(rep("term", length(p.terms)), rep("gene", length(unique(filtered_df$gene_name)))),
                    "is_LA"=c(rep(T, length(p.terms)), #fake it for terms so that conditional fomatting applies only to genes
                              distinct(filtered_df, gene_name, is_LA)$is_LA),
                    "sign_NES"=c(s.nes, rep(0, length(unique(filtered_df$gene_name)))))
  
  ## edges are colored by cluster and thicker if it is a LR
  e.df = filtered_df[,c("term","gene_name","source","cluster","is_LR","is_LA","is_LR.or.LA")]
  
  ## make graph
  g <- graph_from_data_frame(select(e.df, from=term, to=gene_name, cluster), directed=F, vertices=n.df)
  set.seed(123) #have to reset seed every time
  
  c1 = as.character(factor(e.df$cluster, levels=names(cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]]),
                           labels=cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]]))
  c2 = as.character(factor(e.df$cluster, levels=names(cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".light")]]),
                           labels=cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".light")]]))
  E(g)$color <- ifelse(e.df$is_LR, c1, c2)
  E(g)$width <- ifelse(e.df$is_LR, 2, 1)
  E(g)$source <- e.df$source

  c3 = as.character(factor(n.df$sign_NES, levels=c(-1,1,0), labels=c("#CFEBF7","#FFC0B5","grey85")))
  V(g)$color <- ifelse(n.df$is_LA, c3, "grey60")
  
  V(g)$label.cex = ifelse(n.df$type=="term", .5, .7)
  V(g)$label.font = ifelse(n.df$type=="term", 1, 3)
  V(g)$label.color = "black"
  V(g)$label = n.df$label
  
  lay1 = layout_with_kk(g)
  set.seed(123)
  
  return(list("nodes"=n.df, "edges"=e.df, "igraph"=g, "layout"=lay1))
}
