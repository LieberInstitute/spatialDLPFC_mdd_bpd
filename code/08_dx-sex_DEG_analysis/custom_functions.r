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


term2GeneIGRAPH <- function(leadingEdge_DF, source=FALSE, node_size=NA, layout_style=c("kk","fr"), text_title=NULL,
                            color_by=c("term","gene")) {
  set.seed(123) 
  
  #check that all in same direction
  st = table(sign(leadingEdge_DF$NES))
  if(dim(st)>1) warning("Not all terms in same direction (mix of depleted and enriched)")
  
  if(!is.na(node_size)) colnames(leadingEdge_DF)[grep(node_size, colnames(leadingEdge_DF))] = "node_size"

  p.terms = tibble::deframe(distinct(leadingEdge_DF, term, name))
  p.terms = sapply(p.terms, function(x) paste(strwrap(x, width=20), collapse="\n"))
  
  
  n.df = data.frame("name"=c(names(p.terms), unique(leadingEdge_DF$gene_name)),
                    "label"=c(p.terms, unique(leadingEdge_DF$gene_name)),
                    "type"=c(rep("term", length(p.terms)), rep("gene", length(unique(leadingEdge_DF$gene_name))))
  )
  if(color_by=="term") {
    s.nes = tibble::deframe(distinct(leadingEdge_DF, term, sign(NES)))
    #make sure that order is the same
    n.df$sign_NES <- c(s.nes[n.df$label[n.df$type=="term"]], 
                       rep(0, length(unique(leadingEdge_DF$gene_name))))
  }
  if(color_by=="gene") {
    s.nes = tibble::deframe(distinct(leadingEdge_DF, gene_name, sign(NES)))
    n.df$sign_NES <- c(rep(0, length(unique(leadingEdge_DF$term))),
                       #make sure that order is the same
                       s.nes[n.df$label[n.df$type=="gene"]])
  }
  if(!is.na(node_size)) {
    s.size = tibble::deframe(distinct(leadingEdge_DF, gene_name, node_size))
    n.df$node_size <- c(rep(6, length(unique(leadingEdge_DF$term))),
                       #make sure that order is the same
                       s.size[n.df$label[n.df$type=="gene"]])
  }
  ## edges are colored by cluster and thicker if it is a LR
  if(source) {
    e.df = leadingEdge_DF[,c("term","gene_name","source")]
  } else {
    e.df = leadingEdge_DF[,c("term","gene_name")]
  }
  if(color_by=="term") e.df$sign_NES = s.nes[e.df$term]
  if(color_by=="gene") e.df$sign_NES = s.nes[e.df$gene_name]
  ## make graph
  g <- graph_from_data_frame(select(e.df, from=term, to=gene_name), directed=F, vertices=n.df)
  set.seed(123) #have to reset seed every time
  
  if(source) E(g)$source <- e.df$source
  E(g)$color = as.character(factor(e.df$sign_NES, levels=c(-1,1,0), labels=c("skyblue","tomato","grey")))
  
  V(g)$color = as.character(factor(n.df$sign_NES, levels=c(-1,1,0), labels=c("#CFEBF7","#FFC0B5","grey85")))
  #V(g)$label.cex = ifelse(n.df$type=="term", .5, .7)
  V(g)$label.font = ifelse(n.df$type=="term", 1, 3)
  V(g)$label.color = "black"
  V(g)$label = n.df$label

  if(!is.na(node_size)) {
    V(g)$size = n.df$node_size
    V(g)$label.cex = ifelse(n.df$type=="term" | n.df$node_size==3, .5, .7)
  } else {
    V(g)$label.cex = ifelse(n.df$type=="term", .5, .7)
  }

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
