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

generateIGRAPH <- function(node_df, edge_df, size_var, text_title=NULL) {
  set.seed(123) # this is needed for consistent layouts
  # edge data frame
  e.df= select(edge_df, from= reference, to= query, jc= jaccard)
  n.df = select(node_df, name=term, label=name, NES, size_sig=any_of(size_var))

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
                                    labels=c("grey","#CFEBF7","grey","#FFC0B5"))
  )
  V(g)$label.cex = .5
  V(g)$label.color = "black"
  V(g)$label = n.df$label
  V(g)$size = n.df$size_scaled
  
  if(is.null(text_title)) {text_title="Use text_title arg to set title."}
  
  outList = list(igraph= g,
                 layout= layout_with_fr(g, weights=E(g)$jc*3),
                 title_text= text_title,
                 sub_text= paste("Node size linearly corresponds to the number of LA DEGs within leading edge.\nSmallest node =", 
                                 min(n.df$size_sig), "DEG(s), largest node =", max(n.df$size_sig),"DEGs"))
  set.seed(123)
  return(outList)
}

term2GeneIGRAPH <- function(leadingEdge_DF, source=FALSE, layout_style=c("kk","fr"), text_title=NULL) {
  set.seed(123) 
  
  #check that all in same direction
  st = table(sign(leadingEdge_DF$NES))
  if(dim(st)>1) warning("Not all terms in same direction (mix of depleted and enriched)")
  
  p.terms = tibble::deframe(distinct(leadingEdge_DF, term, name))
  p.terms = sapply(p.terms, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  s.nes = tibble::deframe(distinct(leadingEdge_DF, term, sign(NES)))
  
  n.df = data.frame("name"=c(names(p.terms), unique(leadingEdge_DF$gene_name)),
                    "label"=c(p.terms, unique(leadingEdge_DF$gene_name)),
                    "type"=c(rep("term", length(p.terms)), rep("gene", length(unique(leadingEdge_DF$gene_name)))),
                    "sign_NES"=c(s.nes, rep(0, length(unique(leadingEdge_DF$gene_name)))))
  
  ## edges are colored by cluster and thicker if it is a LR
  if(source) {
    e.df = leadingEdge_DF[,c("term","gene_name","source")]
  } else {
    e.df = leadingEdge_DF[,c("term","gene_name")]
  }
  
  ## make graph
  g <- graph_from_data_frame(select(e.df, from=term, to=gene_name), directed=F, vertices=n.df)
  set.seed(123) #have to reset seed every time
  
  if(source) E(g)$source <- e.df$source
  
  V(g)$color = as.character(factor(n.df$sign_NES, levels=c(-1,1,0), labels=c("#CFEBF7","#FFC0B5","grey85")))
  V(g)$label.cex = ifelse(n.df$type=="term", .5, .7)
  V(g)$label.font = ifelse(n.df$type=="term", 1, 3)
  V(g)$label.color = "black"
  V(g)$label = n.df$label
  
  if(layout_style=="kk") {
    lay1 = layout_with_kk(g)
    set.seed(123)
  }
  if(layout_style=="fr") {
    lay1 = layout_with_fr(g)
    set.seed(123)
  }
  
  if(is.null(text_title)) {text_title="Use text_title arg to set title."}

  return(list("igraph"=g, "layout"=lay1, "title_text"=text_title))
}
