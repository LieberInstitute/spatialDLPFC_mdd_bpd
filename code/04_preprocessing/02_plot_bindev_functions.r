plotBiasedFeatures <- function(dataframe, n.features) {
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
