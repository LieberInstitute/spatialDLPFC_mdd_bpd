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

plotCorHeatmap <- function(query_list_level, query_stats, 
                           reference_list_level, reference_stats,
			   .coord_flip = FALSE) {
  both.genes = intersect(top50.both.list[[query_list_level]], rownames(reference_stats))
  cor.res = cor(query_stats[both.genes,],reference_stats[both.genes,])
  cor.df = tibble::rownames_to_column(as.data.frame(cor.res), var="clusters") %>%
    tidyr::pivot_longer(all_of(colnames(cor.res)), names_to=reference_list_level, 
                        values_to="pearson_r")
  cor.df$reference = factor(cor.df[[reference_list_level]], levels= rev(orderList[[reference_list_level]]))
  cor.df$clusters = factor(cor.df$clusters, levels= orderList[[query_list_level]])
  #if want to coord flip also have to change labels so that r_labels aren't flipped and q_labels are
  if(orderList[[query_list_level]][[1]]=="Micro.Vasc") {
    q_labels = c("M/V", orderList[[query_list_level]][-1])
    if(.coord_flip) q_labels = rev(q_labels)
  } else {
      q_labels = orderList[[query_list_level]]
      if(.coord_flip) q_labels = rev(q_labels)
    }
  if(orderList[[reference_list_level]][[1]]=="Micro.Vasc") {
    r_labels = c("M/V", orderList[[reference_list_level]][-1])
    if(!.coord_flip) r_labels = rev(r_labels)
  } else {
      r_labels = orderList[[reference_list_level]]
      if(!.coord_flip) r_labels = rev(r_labels)
    }
  tmp <- ggplot(cor.df, aes(x=clusters, y=reference, fill=pearson_r))+
    geom_tile(color="grey50", linewidth=.3)+
    scale_fill_gradientn("Pearson\nrho", 
                         colors=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(20),
                         limits=c(-max(cor.res), max(cor.res)))+
    scale_x_discrete(query_list_level, labels=q_labels)+
    scale_y_discrete(reference_list_level, labels=r_labels)+
    labs(title=paste("n=",length(both.genes),"genes"),
         subtitle=paste("(of", length(top50.both.list[[query_list_level]]), "top", query_list_level,"genes)"))+
    theme_minimal()+theme(panel.grid= element_blank(), aspect.ratio=1,
                          legend.key.size= unit(10, "pt"), legend.title = element_text(size=8),
                          legend.text = element_text(size=7),
                          plot.title=element_text(size=10), plot.subtitle = element_text(size=8),
                          axis.text.x= element_text(angle=45, hjust=1))
  if(.coord_flip) return(tmp+coord_flip())
  return(tmp)
}

splitDotPlot <- function(df_list_level, df_list, .fill_palette, .fill_name, .y_title, .grob_prop) {
  #select dframe
  tmp = df_list[[df_list_level]]
  #make sure that size legend title reflects true observations
  obs = "spots"
  if(df_list_level=="SZBDMulti-seq") obs="nuclei"
  #swap out Micro.Vasc for M/V in x labels
  x_labs = levels(tmp$clusters)
  if(df_list_level!="PRECAST (smoothed)") x_labs[[1]] = "M/V"
  
  p1 <- ggplot(filter(tmp, gene_type=="protein_coding"), 
               aes(x=clusters, y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
    geom_tile(data=filter(tmp, gene_type=="protein_coding"), aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
    scale_fill_manual(paste0(.fill_name, "\ntop 10\ngenes"), values=.fill_palette)+
    geom_count()+scale_color_gradient("Avg. expr.\n(scaled)", low="white", high="black")+
    scale_size(paste0("Prop. of\n",obs), range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
    scale_x_discrete(labels=x_labs)+
    labs(y=.y_title, subtitle="protein coding genes",
         title=df_list_level)+
    theme_minimal()+theme(axis.text.y=element_text(face="italic"), axis.title.x=element_blank(),
                          legend.title= element_text(size=9),
                          legend.key.size= unit(10, "pt"), plot.title.position = "plot",
                          axis.title.y=element_text(margin=margin(0,20,0,0,"pt")))
  #return(p1)
  p2 <- ggplot(filter(tmp, gene_type!="protein_coding"), 
               aes(x=clusters, y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
    geom_tile(data=filter(tmp, gene_type!="protein_coding"), aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
    scale_fill_manual("", values=.fill_palette)+
    geom_count()+scale_color_gradient("Avg. expr.\n(scaled)", low="white", high="black", guide="none")+
    scale_size(paste0("Prop. of\n",obs), range=c(1,6), limits=c(0,1), breaks=c(0,.5,1), guide="none")+
    scale_x_discrete(labels=x_labs)+
    labs(y="", subtitle="NOT protein coding genes")+
    theme_minimal()+theme(axis.text.y=element_text(face="italic"), axis.title.x=element_blank(),
                          legend.key.size= unit(10, "pt"),
                          axis.title.y=element_text(margin=margin(0,20,0,0,"pt")))
  prop1 = rep(1, length.out=.grob_prop*10)
  prop2 = rep(2, length.out = 10-length(prop1))
  laymat = matrix(c(prop1, prop2))
  return(arrangeGrob(p1, p2, layout_matrix=laymat))
}

dotplotList <- function(top_marker_source, grob_prop) {
  y_title = paste(top_marker_source, "top markers")
  list1 = top.df.list[[top_marker_source]]
  #set fill palette
  if(top_marker_source=="SZBDMulti-seq") fill_palette = cpList$low.res.light
  if(top_marker_source=="PRECAST (smoothed)") fill_palette = cpList$smoothed.light
  if(top_marker_source=="MBv label transfer") fill_palette = cpList$transfer.light
  grobList <- lapply(names(list1), function(x) splitDotPlot(x, list1, fill_palette, 
                                                            .fill_name= gsub(" |-","\n",top_marker_source),
                                                            y_title, grob_prop))
  return(grobList)
}
