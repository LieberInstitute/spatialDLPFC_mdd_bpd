### takes input data for 1 term, 1 dx*sex group, 1 LR cluster and outputs dataframes with attr for plotting manhattan and manhattan-labeled volcanoes
formatManhattan <- function(filtered_LR_df, specific_fgsea_results, .gmtList, term) {
  #filtered_LR_df is the compiled results LR DF filtered to 1 dx*sex group and 1 cluster
  #specific_fgsea_results is the list level of the fgsea results list for 1 dx*group and 1 cluster
  lf.df = group_by(filtered_LR_df, gene_name) %>% summarise(avg_lf=mean(logFC)) %>% arrange(desc(avg_lf))
  lf.order = lf.df$avg_lf
  names(lf.order) = lf.df$gene_name
  
  ledge = tibble::deframe(specific_fgsea_results[,c("pathway","leadingEdge")])
  
  plot.df = data.frame("term"=rep(term, length(lf.order)), 
                       "gene_name"=names(lf.order),
                       "in_term"=names(lf.order) %in% .gmtList[[term]],
                       "in_leadingEdge"=names(lf.order) %in% ledge[[term]],
                       "logFC"=lf.order,
                       "rank"=1:length(lf.order), row.names = NULL)
  
  limitsList = list("max_rank"=max(plot.df$rank),
                    "max_logFC"=max(plot.df$logFC),
                    "min_logFC"=min(plot.df$logFC))
  plot.df = filter(plot.df, in_term==T)
  attr(plot.df, 'limitsList') = limitsList
  
  term_results = specific_fgsea_results[specific_fgsea_results$pathway==term,]
  labelList = list("term"=term,
                   "NES"=term_results[["NES"]], 
                   "n_leadingEdge"=length(ledge[[term]]),
                   "n_geneSet_overlap"= term_results[["size"]],
                   "n_geneSet_total"=length(.gmtList[[term]]),
                   "padj"=term_results[["padj"]])
  labelList$label = paste0("NES= ", round(labelList$NES, 2), "\n(", 
                           labelList$n_leadingEdge, "/ ", labelList$n_geneSet_overlap, ")\n",
                           "padj= ", format(labelList$padj, scientific=T, digits=2))
  attr(plot.df, 'labelList') = labelList
  
  plot.df2 = mutate(filtered_LR_df, in_term = gene_name %in% .gmtList[[term]],
                    in_leadingEdge = gene_name %in% ledge[[term]],
                    point_col= factor(paste(in_term, in_leadingEdge), 
                                      levels=c("FALSE FALSE", "TRUE FALSE", "TRUE TRUE"), 
                                      labels=c("not present", "in gene set\nfor term", "in leading\nedge")))
  attr(plot.df2, 'labelList') = labelList
  
  return(list("manhattan"=plot.df, "volcano"=plot.df2))
}


### returns ggplot object of manhattan plot
plotManhattan <- function(man_df) {
  p1 <- ggplot(man_df, aes(x=rank, xend=rank, y=0, yend=logFC, color=in_leadingEdge))+
    geom_segment(linewidth=.3)+scale_color_manual(values=c("grey","black"), guide="none")+
    scale_x_continuous(limits=c(0, attr(man_df, 'limitsList')$max_rank), expand=c(.01,.01))+
    geom_text(aes(x=median(c(0, attr(man_df, 'limitsList')$max_rank)), y=.1, 
                  label=attr(man_df, "labelList")$label), 
              vjust=0, size=3, color="black")+
    theme_minimal()+labs(title=attr(man_df, "labelList")$term, y="logFC")+
    theme(axis.title=element_text(size=8), axis.text=element_text(size=7),
          plot.title=element_text(size=10, face="bold"))
  return(p1)
}


### returns ggplot object of volcano plot labeled with manhattan leading edge genes
plotVolcano <- function(vol_df) {
  v1 <- ggplot(vol_df, aes(x=logFC, y=-log10(adj.P.Val), color=point_col))+
    geom_point(data=filter(vol_df, in_term==F), size=.5)+
    geom_point(data=filter(vol_df, in_term==T), size=.5)+
    geom_hline(aes(yintercept=-log10(.05)), lty=2, color="red3")+
    ggrepel::geom_label_repel(data=filter(vol_df, in_term==T, adj.P.Val<.05),
                              aes(label=gene_name), min.segment.length = 0, size=3, show.legend = FALSE)+
    geom_text(data=data.frame("x1"=round(max(abs(vol_df$logFC))/2, 1), "y1"=-log10(.05), "lab1"="adj. p<.05"),
              aes(x=x1, y=y1, label=lab1), color="red3", size=3, fontface="italic", 
              vjust=0, hjust=0, nudge_y = .05)+
    scale_color_manual("Gene status:", values=c("grey80","grey50","black"),
                       guide=guide_legend(override.aes = list(size=2)))+
    xlim(-max(abs(vol_df$logFC)), max(abs(vol_df$logFC)))+
    labs(title=attr(vol_df, 'labelList')$term,
         subtitle=paste(unique(vol_df$group), unique(vol_df$sex), unique(vol_df$cluster)))+
    theme_bw()+theme(text=element_text(size=8), legend.box.spacing= unit(1,"pt"))
  return(v1)
}
