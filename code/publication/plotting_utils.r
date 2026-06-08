suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

# pre-reqs for dotplot, consensus, and violin
source("code/09_DEG_GRN/load_DEGs.r")
all_clusters = c("L-A sm","L-A se","Micro.Vasc se","Astro se","L1 sm",
                 "L2 sm","L2.3 se","L3.4 sm","L4 se",
                 "Inhb se","L5 sm","L5 se","L6 sm","L6 se","WM sm","Oligo se")

# pre-reqs for prop detected boxplot
load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")

# pre-reqs for mean ratio bar plot
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")
sn.col.pal = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
#            "Micro"=cpList$low.res.light[["L3"]],
	    "Micro"="#C28658",
            cpList$low.res.light[c("Astro","Oligo")],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]]
)

# pre-reqs for violin plot
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
col.pal_fill = c("NTC"="#CBCBCB", "MDD"="#FDCFBF", "BPD"="#D3BFE0",
                 "NTC False"="transparent", "MDD False"="transparent", "BPD False"="transparent",
                 "NTC True"=cpList$dx.pal[["NTC"]], "MDD True"=cpList$dx.pal[["MDD"]], "BPD True"=cpList$dx.pal[["BPD"]])

col.pal_color = c("NTC False"="grey50", "MDD False"="grey50", "BPD False"="grey50",
                  "NTC True"="black", "MDD True"="black", "BPD True"="black")


# functions
getConsensus <- function(module_genes, sig.df, return_DF=FALSE) {
	con.df = filter(sig.df, gene_name %in% module_genes)
	con.df$x_labels = factor(con.df$cluster_source, levels=all_clusters,
                         labels=c("L-A","L-A","M.V","Ast","L1",
                                  #"L2","L2.3","L3.4","L4",
                                  "L2/3","L2/3","L3/4","L3/4",
                                  "Inb",
                                  "L5","L5","L6","L6","WM/O","WM/O"))
	#number of consensus to filter to depends on if domain/model result is repeated
	con.df$max_consensus = as.numeric(as.character(factor(con.df$cluster_source, levels=all_clusters,
                                labels=c(2,2,1,1,1,
                                         2,2,2,2,
                                         1,
                                         2,2,2,2,2,2))))
	consensus_genes = group_by(con.df, gene_name, sex.group, dir, x_labels, max_consensus) %>%
	  tally(name="n_consensus") %>% filter(n_consensus==max_consensus) 
	if(return_DF) {
		return(consensus_genes)
	} else {
		return(unique(consensus_genes$gene_name))
	}
}


getDotplot <- function(ordered_genes, de.df, color_scale_limits=c(-2.5,2.5)) {
	dot.df = filter(de.df, gene_name %in% ordered_genes) %>% mutate(source=factor(source, levels=c("sm","se"))) 
	dot.df$is_sig = dot.df$adj.P.Val2<.05
	dot.df$gene_name = factor(dot.df$gene_name, levels=rev(ordered_genes))
	dot.df$x_labels = factor(dot.df$cluster_source, levels=all_clusters,
                         labels=c("L-A","L-A","M.V","Ast","L1",
                                  #"L2","L2.3","L3.4","L4",
                                  "L2/3","L2/3","L3/4","L3/4",
                                  "Inb",
                                  "L5","L5","L6","L6","WM/O","WM/O"))
	dot.df$size2 = as.numeric(as.character(factor(paste(dot.df$source, dot.df$is_sig), 
                                              levels=c("sm FALSE","se FALSE","sm TRUE","se TRUE"),
                                              labels=c(1,1,4,3))))
	p1 <- ggplot(dot.df, aes(x=x_labels, y=gene_name, fill=logFC, size=size2))+
	  geom_count(aes(shape=source, color=is_sig))+
	  scale_shape_manual(values=c("sm"=23, "se"=21))+
	  scale_size_identity("adj. p")+ 
	  scale_color_manual(values=c("FALSE"="grey", "TRUE"="black"))+
	  scale_fill_gradientn("logFC",colors=RColorBrewer::brewer.pal(n=5,"RdBu")[5:1],
	                       limits=color_scale_limits)+
	  facet_grid(cols=vars(sex.group), #rows=vars(gene_group), 
	             scales="free_y", space="free_y")+
	  guides(fill=guide_colorbar(theme=theme(legend.key.height=unit(12,"pt"), legend.key.width=unit(36,"pt"))))+
	  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=8), axis.title.x=element_blank(),
	                   strip.background.y = element_blank(), strip.text.y = element_blank(),
	                   legend.position="bottom", axis.title.y=element_blank(),
	                   legend.text = element_text(size=8))
	return(p1)
}


getDetectedBoxplot <- function(ordered_genes, spe_summ) {
	gids = rownames(spe_summ)[rowData(spe_summ)$gene_name %in% ordered_genes]
	names(gids) = rowData(spe_summ)[gids,"gene_name"]
	stopifnot(length(ordered_genes)==length(gids))

	df1 = as.data.frame(assay(spe_summ, "logcounts.prop.detected")[gids,])
	df2 = tidyr::pivot_longer(tibble::rownames_to_column(df1, var="gene_id"), all_of(colnames(df1)), names_to="sample", values_to="prop.spots.detected") %>%
	  left_join(as.data.frame(rowData(spe_summ)[gids,c("gene_id","gene_name")])) %>%
	  mutate(gene_name= factor(gene_name, levels=rev(ordered_genes)))

	p1.1 <- ggplot(df2, aes(y=gene_name, x=prop.spots.detected))+
	  geom_boxplot(outlier.size=.5, fill="grey")+xlim(0,1)+
	  labs(title=" ", x="prop. spots")+
	  theme_minimal()+theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.text.x=element_text(angle=60, hjust=1, size=8),
                        plot.margin = margin(.2,.2,1.5,.2,"cm"))
	return(p1.1)
}


getMeanRatioBar <- function(ordered_genes, sce_summ) {
	gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% ordered_genes]
	names(gids) = rowData(sce_summ)[gids,"gene_name"]	

	df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
	colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

	bar.df = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), 
                             names_to="cellType", values_to="mean.expr") %>%
	  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")]))

	if(length(ordered_genes)!=length(gids)) {
		add.genes = setdiff(ordered_genes, names(gids))
		for(i in add.genes) {
		  bar.df = add_row(bar.df, gene_id=i, cellType = c("Vasc","Astro","Oligo","Micro","InhN","ExcN"), mean.expr=NA, gene_name=i)
		}
	}

	bar.df$cellType = factor(bar.df$cellType, levels=c("Vasc","Astro","Oligo","Micro","InhN","ExcN"))
	bar.df$gene_name = factor(bar.df$gene_name, levels=rev(ordered_genes))

	p1.2 = ggplot(bar.df, aes(y=gene_name, x=mean.expr, fill=cellType))+
	  geom_bar(stat="identity", position="fill")+
	  scale_fill_manual(values=sn.col.pal)+scale_y_discrete(position="right")+
	  labs(title=" ", x="mean\nexpr")+
	  theme_minimal()+theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
	                        axis.title.y=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
	                        legend.position="bottom", legend.title=element_blank(),
	                        legend.key.size = unit(10,"pt"))
	return(p1.2)
}

getViolin <- function(.plot.genes, .y.upper.bound=11, .y.breaks=c(0,2,4,6,8,10), version=c("seurat","precast")) {
	if(version=="seurat") {
		load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
		return(getViolin_seurat(.plot.genes, .y.upper.bound, .y.breaks, spe_pseudo))
	}
	if(version=="precast") {
		load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
		return(getViolin_precast(.plot.genes, .y.upper.bound, .y.breaks, spe_pseudo))
	}
	if(!version %in% c("seurat","precast")) stop("Please specify pseudobulk version to use as one of 'seurat' or 'precast'.")
	if(length(version)>1) stop("Please specify pseudobulk version to use as *one* of 'seurat' or 'precast'.")
}

getViolin_seurat <- function(plot.genes, y.upper.bound=11, y.breaks=c(0,2,4,6,8,10), spe_pseudo) {
  for (j in plot.genes) {
    colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
  }
  
  summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
    tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=plot.genes),
           seurat_label="L-A")
  summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
    tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=plot.genes))
  all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright)[c(1:4,8,5:7)])))
  
  all.df_filt = filter(all.df, key_genes %in% plot.genes)
  all.df_filt2 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% summarise(ypos=mean(logcounts), ysd=sd(logcounts), n=n(), yse=ysd/sqrt(n)) #%>%
  
  
  
  sig.df = do.call(rbind, sigList[c("se.la","se.lr")]) %>% 
    filter(gene_name %in% plot.genes, group!="MDD.BPD") %>%
    mutate(is_sig=T) %>% tidyr::separate_rows(group, sep="\\.") %>%
    select(gene_name, sex, condition=group, cluster, is_sig)
  
  all.df_filt = left_join(all.df_filt, sig.df, by=c("key_genes"="gene_name","sex","condition","cluster")) %>%
    mutate(is_sig= ifelse(is.na(is_sig), "False", "True"),
           point_color= paste(condition, is_sig),
           condition=factor(condition, levels=c("NTC","MDD","BPD")),
           sex=factor(sex, levels=c("F","M")),
           key_genes=factor(key_genes, levels=plot.genes),
           cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright)[c(1:4,8,5:7)])))
  
  all.df_filt2 = left_join(all.df_filt2, sig.df, by=c("key_genes"="gene_name","sex","condition","cluster")) %>%
    mutate(is_sig= ifelse(is.na(is_sig), "False", "True"),
           point_color= paste(condition, is_sig),
           condition=factor(condition, levels=c("NTC","MDD","BPD")),
           sex=factor(sex, levels=c("F","M")),
           key_genes=factor(key_genes, levels=plot.genes),
           cluster=factor(cluster, levels=c("L-A", names(cpList$transfer.bright)[c(1:4,8,5:7)])))
  
  
  
  p1 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
    geom_violin(aes(fill=condition), scale="width", trim=F, bounds=c(0,y.upper.bound), position=position_dodge(width=.8), 
                color="transparent")+
    facet_grid(cols=vars(sex), rows=vars(key_genes), switch="y",
               labeller= as_labeller(c("F"="Female","M"="Male", plot.genes)))+
    geom_text(data=filter(all.df_filt2, cluster=="L4"), aes(x=cluster, y=1, label=key_genes),
              color="grey50", size=2, fontface="italic", hjust=.5)+
    geom_crossbar(data=all.df_filt2, aes(group=condition, y=ypos, ymax=ypos+ysd, ymin=ypos-ysd, 
                                         fill=point_color, color=point_color),
                  position = position_dodge(width=.8), width=.6, linewidth=.3)+
    scale_color_manual(values=col.pal_color, guide="none")+
    scale_fill_manual(values=col.pal_fill, guide="none")+
    scale_y_continuous(breaks=y.breaks)+coord_cartesian(ylim=c(0,y.upper.bound))+
    scale_x_discrete(labels=c("L-A","M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))+
    theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                     text=element_text(size=8), axis.text=element_text(size=6),
                     axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks = element_line(linewidth=.2), strip.text.y.left = element_blank(),
                     panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2))
  return(p1)
}

getViolin_precast <- function(plot.genes, y.upper.bound=11, y.breaks=c(0,2,4,6,8,10), spe_pseudo) {
  for (j in plot.genes) {
    colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
  }
  
  summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
    tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=plot.genes),
           smoothed_k9_1663="L-A")
  summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "smoothed_k9_1663", plot.genes)]) %>%
    tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
    mutate(key_genes= factor(key_genes, levels=plot.genes))
  all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(smoothed_k9_1663, levels=c("L-A", names(cpList$smoothed.bright))))
  
  all.df_filt = filter(all.df, key_genes %in% plot.genes)
  all.df_filt2 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% summarise(ypos=mean(logcounts), ysd=sd(logcounts), n=n(), yse=ysd/sqrt(n)) #%>%
  
  
  
  sig.df = do.call(rbind, sigList[c("sm.la","sm.lr")]) %>% 
    filter(gene_name %in% plot.genes, group!="MDD.BPD") %>%
    mutate(is_sig=T) %>% tidyr::separate_rows(group, sep="\\.") %>%
    select(gene_name, sex, condition=group, cluster, is_sig)
  
  all.df_filt = left_join(all.df_filt, sig.df, by=c("key_genes"="gene_name","sex","condition","cluster")) %>%
    mutate(is_sig= ifelse(is.na(is_sig), "False", "True"),
           point_color= paste(condition, is_sig),
           condition=factor(condition, levels=c("NTC","MDD","BPD")),
           sex=factor(sex, levels=c("F","M")),
           key_genes=factor(key_genes, levels=plot.genes),
           cluster=factor(smoothed_k9_1663, levels=c("L-A", names(cpList$smoothed.bright))))
  
  all.df_filt2 = left_join(all.df_filt2, sig.df, by=c("key_genes"="gene_name","sex","condition","cluster")) %>%
    mutate(is_sig= ifelse(is.na(is_sig), "False", "True"),
           point_color= paste(condition, is_sig),
           condition=factor(condition, levels=c("NTC","MDD","BPD")),
           sex=factor(sex, levels=c("F","M")),
           key_genes=factor(key_genes, levels=plot.genes),
           cluster=factor(cluster, levels=c("L-A", names(cpList$smoothed.bright))))
  
  
  
  p1 <- ggplot(all.df_filt, aes(x=cluster, y=logcounts))+
    geom_violin(aes(fill=condition), scale="width", trim=F, bounds=c(0,y.upper.bound), position=position_dodge(width=.8), 
                color="transparent")+
    facet_grid(cols=vars(sex), rows=vars(key_genes), switch="y",
               labeller= as_labeller(c("F"="Female","M"="Male", plot.genes)))+
    geom_text(data=filter(all.df_filt2, cluster=="L3.4"), aes(x=cluster, y=1, label=key_genes),
              color="grey50", size=2, fontface="italic", hjust=.5)+
    geom_crossbar(data=all.df_filt2, aes(group=condition, y=ypos, ymax=ypos+ysd, ymin=ypos-ysd, 
                                         fill=point_color, color=point_color),
                  position = position_dodge(width=.8), width=.6, linewidth=.3)+
    scale_color_manual(values=col.pal_color, guide="none")+
    scale_fill_manual(values=col.pal_fill, guide="none")+
    scale_y_continuous(breaks=y.breaks)+coord_cartesian(ylim=c(0,y.upper.bound))+
    scale_x_discrete(labels=c("L-A","L1","L2","L3.4","L5","L6","WM"))+
    theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                     text=element_text(size=8), axis.text=element_text(size=6),
                     axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks = element_line(linewidth=.2), strip.text.y.left = element_blank(),
                     panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2))
  return(p1)
}
