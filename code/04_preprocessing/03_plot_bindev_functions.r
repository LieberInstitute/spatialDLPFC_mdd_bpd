plotBinDevResults <- function(dataframe, n.features) {
	mean1 = mean(dataframe$r.diff)
	sd.outlier = sd(dataframe$r.diff)
	m1 = ceiling(max(dataframe$r.diff)/(mean1+sd.outlier))
	sd.list=c(1:(m1+1))*sd.outlier
	break.list = mean1+sd.list

	dataframe$all.slide.outlier = cut(abs(dataframe$r.diff), breaks=c(0,break.list),labels=FALSE, include.lowest=TRUE)*sign(dataframe$r.diff)

	dataframe$all.slide.outlier.color = cut(abs(dataframe$all.slide.outlier), breaks=seq(0,max(dataframe$all.slide.outlier)+5, by=5), include.lowest=TRUE)

	colpal = c("grey90",viridis::viridis_pal(option="turbo")(length(levels(dataframe$all.slide.outlier.color))-1))

	# SCATTER (all)
	p1 <- ggplot(dataframe, aes(x=rank_default, y=rank_brain, color=all.slide.outlier.color))+
		geom_point()+
		scale_color_manual(values=colpal)+
		labs(title=paste("Binomial deviance;",n.features,"features"), subtitle="color based on SD of pooled slides", x="rank_default", y="rank_brain")+
		theme_bw()+theme(legend.position="none")

	# HISTOGRAM
	p2 <- ggplot(dataframe, aes(all.slide.outlier, fill=all.slide.outlier.color))+
		geom_histogram(color="grey")+
		scale_fill_manual(values=colpal)+
		scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
			breaks=10^(0:ceiling(log10(max(table(dataframe$all.slide.outlier))))))+
		labs(title=paste("Binomial deviance;",n.features,"features"), subtitle="SD(rank diff) [(rank_brain)-(rank_default)]", 
			x="n SD from mean rank diff", y="# spots")+
		theme_bw()+theme(legend.position="none")

	# SCATTER (by slide)
	p3 <- ggplot(dataframe, aes(x=rank_default, y=rank_brain))+
		geom_point(data=filter(dataframe, all.slide.outlier<=5), color="grey90")+
		geom_point(data=filter(dataframe, all.slide.outlier>5), aes(color=all.slide.outlier.color))+
		scale_color_manual(values=colpal[2:length(colpal)])+
		facet_wrap(vars(slide), ncol=2)+
		labs(color="n SD(rank diff)")+
		theme_bw()

	lay=rbind(c(1,3,3),c(2,3,3))
	list(gridExtra::grid.arrange(p1, p2, p3, layout_matrix=lay), dataframe)
}

plotBiasedFeatures <- function(dataframe, sd.cutoff, spatialObject) {
	plot.genes = unique(filter(dataframe, all.slide.outlier>sd.cutoff)$gene)
	names(plot.genes) = rowData(spatialObject)[plot.genes,"gene_name"]

	best.rank.slide.df = filter(dataframe, gene %in% plot.genes) %>% group_by(slide, gene, gene_name) %>% 
		summarize(best.rank.slide=min(rank_brain), .groups="drop")
	best.rank.all.df = filter(dataframe, gene %in% plot.genes) %>% group_by(gene, gene_name) %>% 
		summarize(best.rank.all = min(rank_brain), .groups="drop") %>%
		arrange(best.rank.all) %>% mutate(index=row_number(), i2=length(plot.genes)-index, ytext=paste(best.rank.all,gene_name, sep=" - ")) %>%
		arrange(i2) %>% mutate(ylabel=factor(i2, levels=i2,labels=ytext))

	plot.rank.df = filter(dataframe, gene %in% plot.genes) %>% select(slide, gene, gene_name, rank_brain, rank_default, r.diff) %>%
		left_join(best.rank.all.df, by=c("gene","gene_name"))

	p1 <-  ggplot(plot.rank.df, aes(x=slide, y=ylabel, size=rank_default, color=rank_brain))+
		geom_count()+scale_color_viridis_c(direction=-1)+
		scale_size(breaks=c(100,500,1000,1500,2000))+
		labs(x="slide",y="(best rank across slides) - gene name", size="rank (no batch)", color="rank (batch = sample)",
			subtitle="missing spot indicates gene not in top deviant genes for default or subject-batch")+
		theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

	l2 = unique(spatialObject$sample_id)
	names(l2) = lapply(l2, function(x) unique(colData(spatialObject)[spatialObject$sample_id==x,"brain"]))
	l2 = lapply(l2, function(x) spatialObject[,colData(spatialObject)$sample_id==x])

	plot.genes.df = do.call(rbind, lapply(seq_along(l2), function(x) {
		m1 = as.data.frame(do.call(rbind,
			lapply(plot.genes, function(y) {
				gex = logcounts(l2[[x]])[y,]
				c("avg"=mean(gex), "n"=sum(gex>0)/length(gex))
			})
		))
		m1$gene_name = names(plot.genes)
		m1$gene = plot.genes
		m1$brain = names(l2)[x]
		m1$slide = unique(l2[[x]]$slide)
		m1$position = unique(l2[[x]]$position)
		m1$condition = unique(l2[[x]]$condition)
		return(m1)
	}))

	plot.genes.df <- left_join(plot.genes.df, best.rank.slide.df, by=c("slide","gene","gene_name")) %>%
		left_join(best.rank.all.df, by=c("gene","gene_name")) %>%
		mutate(xlabel=paste(position, brain, condition))

	p2 <- ggplot(plot.genes.df, aes(x=xlabel, y=ylabel, size=n, color=avg))+
		geom_count()+scale_color_viridis_c(option="F", direction=-1)+
		facet_wrap(vars(slide), ncol=6, scales="free_x")+
		labs(x="", y="(best rank across slides) - gene name",size="expr (prop. spots)", color="expr (avg. logcounts)")+
		theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1))

	return(list(p1, p2))
}
