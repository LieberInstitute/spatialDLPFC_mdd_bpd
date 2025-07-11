dotplotDF <- function(source_sce, gene_id_list, summarize_groups=FALSE, cluster_labels=NULL) {
  dotplot.df = left_join(tibble::rownames_to_column(as.data.frame(assay(source_sce, "logcounts.mean")[gene_id_list,]), var="gene_id") %>%
                           tidyr::pivot_longer(colnames(source_sce), names_to="clusters", values_to="mean_expr"),
                         tibble::rownames_to_column(as.data.frame(t(scale(t(assay(source_sce, "logcounts.mean")[gene_id_list,])))), var="gene_id") %>%
                           tidyr::pivot_longer(colnames(source_sce), names_to="clusters", values_to="mean_expr_scaled")) %>%
    left_join(tibble::rownames_to_column(as.data.frame(assay(source_sce, "logcounts.prop.detected")[gene_id_list,]), var="gene_id") %>%
                tidyr::pivot_longer(colnames(source_sce), names_to="clusters", values_to="prop_spots")) %>%
    left_join(rdata[,c("gene_id","gene_name","gene_type")], by="gene_id")
  if(!summarize_groups) {return(dotplot.df)}
  else {
    if(is.null(cluster_labels)) {error("To summarize by groups, name of colData column of clusters to summarize to must be provided as 'cluster_labels'.")}
    cdata = as.data.frame(colData(source_sce))
    cdata$new_clusters = cdata[[cluster_labels]]
    dotplot.df2 = left_join(dotplot.df, cdata[,c("sample_id","new_clusters")],
                            by=c("clusters"="sample_id")) %>%
      group_by(new_clusters, gene_name, gene_id, gene_type) %>% summarise_at(c("mean_expr_scaled","prop_spots"), mean)
    colnames(dotplot.df2)[grep("new_clusters", colnames(dotplot.df2))] = "clusters"
    return(dotplot.df2)
  }
}

getTopGenes <- function(enrich_stats, top_n=100, return_matrix=FALSE) {
  rank.mtx = apply(enrich_stats, 2, rank)
  rank.mtx = (nrow(rank.mtx)+1)-rank.mtx
  top.ranks = rank.mtx<=top_n
  #warn if top genes include down-regulated genes
  if(sum(enrich_stats[top.ranks]<0)>0) warning(paste("top_n =", top_n, "produces a gene set where some top genes have negative test statistics."))
  #return rownames or matrix if requested (matrix allows to easily map top genes to cluster)
  if(return_matrix) {
    return(top.ranks)
  } else {
    top.genes = rownames(top.ranks)[rowSums(top.ranks)>0]
    return(unique(top.genes))
  }
}


VolcanoPlot <- function(source_df, cluster_labels) {
	#, gene_id_list) {
  y_max = ceiling(max(-log10(source_df$adj.P.Val)))
  if(y_max==Inf) {y_max = ceiling(max(-log10(source_df$adj.P.Val[source_df$adj.P.Val!=0])))}
  y_limits = c(0,y_max)
  x_limits = c(-ceiling(max(abs(source_df$logFC))),ceiling(max(abs(source_df$logFC))))
  
  clus = levels(source_df[[cluster_labels]])
  names(clus) = clus
  plist = lapply(clus, function(x) {
    tmp = source_df[source_df[[cluster_labels]]==x,]
    ggplot(tmp, aes(x=logFC, y=-log10(adj.P.Val), color=gene_type_ptn))+
      geom_point(size=.3, color="grey")+#scale_color_manual(values=c("grey","black"), guide="none")+
      #geom_point(data=filter(tmp, logFC>0) %>% slice_min(n=10, adj.P.Val), size=.5, color="black")+
      geom_point(data=slice_max(tmp, n=50, t), size=.5, color="black")+
      geom_point(data=slice_max(tmp, n=50, logFC), size=.5, color="black")+
      #geom_point(data=filter(tmp, gene_id %in% gene_id_list), size=.5, color="black")+
      xlim(x_limits)+#ylim(y_limits)+
      labs(title=x)+
      theme_bw()+theme(panel.grid.minor=element_blank(), plot.title=element_text(size=9),
                       axis.title.x=element_blank(), axis.title.y=element_blank())
  })
  return(plist)
}
