lossPlot = function(gene_set, k_clusters, jobid=NA) {
  #gene_set = character specifying gene set (starting with 'n' for most instances except for H-M gene set)
  #k_clusters = numeric
  #jobid = optional job specification
  x = grep(paste0(".*",gene_set), list.files("code/05_clustering/PRECAST/logs"), value=T)
  x = grep(paste0(".*k", k_clusters), x, value=T)
  #remove plotting logs
  x = x[grep("^precast_", x)]
  #quit if multiple files
  if(!is.na(jobid)) {x= grep(paste0(".*",jobid), x, value=T)}
  if(length(x)>1) {
    if(is.na(jobid)) {stop(paste0("More than 1 log file found for gene set ", gene_set," and k=", k_clusters,". Set jobid and retry."))}
    if(!is.na(jobid)) {stop(paste0("Log for gene_set ", gene_set, " and k=", k_clusters, " with jobid=", jobid, " not found. Check for typos and try again."))}
    }
  cat("Plotting for file:", x,"\n")
  x = readLines(paste0("code/05_clustering/PRECAST/logs/", x))
  
  iters = x[grep("iter =", x)]
  iter.list = strsplit(iters, ", ")
  df = data.frame(iter=2:(length(iters)+1), 
                  loglik=as.numeric(sapply(iter.list, function(x) substr(x[[2]], start=9, stop=20))),
                  d.loglik=as.numeric(sapply(iter.list, function(x) substr(x[[3]], start=9, stop=16))))
  
  #modified to be saved with annotation heatmap
  #png(file=paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters, "_lglk-loss-plot.png"), bg="white")
  plot(df[["iter"]], -df[["loglik"]], main=paste0("PRECAST ", gene_set, " k=", k_clusters), 
       xlab="iteration", ylab="neg. loglik")
  #dev.off()
  #cat("\nLoss plot saved to:", paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters, "_lglk-loss-plot.png"),"\n")
}

updateColData <- function(object, gene_set, k_clusters) {
  #object = spe
  #gene_set = character specifying gene set (starting with 'n' for most instances except for H-M gene set)
  #k_clusters = numeric
  if(substr(gene_set, start=0, stop=1)=="n") {gene_set2 = substr(gene_set, start=2, stop=6)} 
  if(gene_set=="H-M-markers") {gene_set2 = "HM"}
  
  load(paste0("processed-data/05_clustering/PRECAST/srt_precast_k-", k_clusters, "_", gene_set,".Rdata"))
  
  object$seurat_key = paste(object$sample_id, colnames(object), sep="_")
  mdata = seuInt@meta.data[object$seurat_key,]
  stopifnot(identical(rownames(mdata), object$seurat_key))
  colData(object)[[paste0("precast_k", k_clusters,"_", gene_set2)]] = as.factor(mdata$cluster)
  
  #quickResaveHDF5SummarizedExperiment(object)
  #cat("\nPRECAST clusters", gene_set, "genes, k=", k_clusters, "updated to spe with quickResave\n")
  return(object)
}

generateSpotPlots <- function(object, gene_set, k_clusters, color_palette) {
  #object = spe
  #gene_set = character specifying gene set (starting with 'n' for most instances except for H-M gene set)
  if(substr(gene_set, start=0, stop=1)=="n") {gene_set2 = substr(gene_set, start=2, stop=6)} 
  if(gene_set=="H-M-markers") {gene_set2 = "HM"}
  cluster_col = paste0("precast_k", k_clusters,"_", gene_set2)
  #k_clusters = numeric
  #color_palette = named vector of colors with names being CLUSTER NUMBER

  object$dummy_slide = ifelse(object$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", object$slide)
  object$array2 = ifelse(object$slide=="V13B23-283", "A1", object$array)
  
  seed = levels(as.factor(object$dummy_slide))
  slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
  slideList = lapply(slideList, function(x) {
    do.call(cbind, lapply(x, function(y) object[,object$dummy_slide==y]))
  })
  
  clusPlot = lapply(slideList, function(x) {
    suppressMessages(
      plotSpots(x, annotate=cluster_col, 
                point_size=.5, sample_id="sample_id")+
        scale_color_manual(values=color_palette)+
        facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
        theme(panel.background=element_rect(fill="grey30"))
    )
  })
  
  return(clusPlot)
  #pdf(file=paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters,"_spot-plots.pdf"), width=12, height=16)
  #clusPlot[[1]]
  #clusPlot[[2]]
  #clusPlot[[3]]
  #clusPlot[[4]]
  #clusPlot[[5]]
  #dev.off()
  #cat("\nSaved spot plots to:",paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters,"_spot-plots.pdf"), "\n")
}

annotationHeatmap <- function(object, gene_set, k_clusters, layer_marker_df) {
  #object = spe
  #gene_set = character specifying gene set (starting with 'n' for most instances except for H-M gene set)
  if(substr(gene_set, start=0, stop=1)=="n") {gene_set2 = substr(gene_set, start=2, stop=6)} 
  if(gene_set=="H-M-markers") {gene_set2 = "HM"}
  cluster_col = paste0("precast_k", k_clusters,"_", gene_set2)
  #k_clusters = numeric
  valid.genes = intersect(rownames(object), layer_marker_df$gene_id)
  cat("\nNumber of valid genes =", length(valid.genes), "\n")
  valid.genes[1:4]
  layer_marker_df = layer_marker_df[layer_marker_df$gene_id %in% valid.genes,]
  dim(layer_marker_df)
  spe_summ = aggregateAcrossCells(object[layer_marker_df$gene_id,], 
                                  #it wouldn't work unless i kep sample_id in coldata so needed to add extra averaging step below
                                  ids=colData(object)[,c("sample_id",cluster_col)], 
                                  statistics=c("mean"),
                                  use.assay.type="logcounts")
  #average together all sample_ids
  prc_clus = 1:k_clusters
  names(prc_clus) = paste0("c",prc_clus)
  m1 = sapply(prc_clus, function(x) {
    idx = colData(spe_summ)[[cluster_col]]==x
    rowMeans(logcounts(spe_summ)[,idx])
  })
  #format layer marker rows
  row_annot = data.frame("cluster"=factor(layer_marker_df$domain_simple, levels=c("GABA","Vasc","L1","L2","L3","L4","L5","L6","WM (1)","WM (2)")))
  rownames(row_annot) = layer_marker_df$gene_id
  annot_colors=list("cluster"=c("Vasc"="#FF7F00","L1"= "#1F78B4", "L2"="#33A02C", "L3"="#A6CEE3", "L4"="grey", "L5"="#B2DF8A", 
	"L6"="#FB9A99", "WM (1)"="#8B0000", "WM (2)"="#E31A1C", "GABA"="#CAB2D6"))
  #plot and save heatmap
  hmp <- pheatmap(m1[rownames(row_annot)[order(row_annot[,1])],], 
                 annotation_row=row_annot, annotation_colors = annot_colors, 
                 show_rownames = F, cluster_rows=F, annotation_names_row = F,
                 scale="row", angle_col = 0, treeheight_col = 10, silent=T)
  #modified to be returned and saved with loss plot
  return(hmp[[4]])
  #ggsave(paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters, "_initial-layer-annotation.png"), hmp[[4]], 
  #       bg="white", width=6, height=7, units="in")
  #cat("\nInitial layer annotation heatmap saved to:",paste0("plots/05_clustering/PRECAST/PRECAST_", gene_set, "-k", k_clusters, "_initial-layer-annotation.png"),"\n")
}
