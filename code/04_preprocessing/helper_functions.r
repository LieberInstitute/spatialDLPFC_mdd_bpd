dotplotDF <- function(plot.genes, spe, norm.to.mbp=TRUE, order.by.rank=FALSE, rank.df=NULL, ...) {
	#description
	### plot.genes: list of ensembl IDs of genes for which expression will be extracted from spe object
	### spe: spatial experiment object with logcounts assay for plotting gene expression
	### norm.to.mbp: logical indicating whether or not to divide image avg gene expression by mbp
	### order.by.rank: logical indicating if the y axis of the resulting dotplot should be organized by genes with decreasing max rank (if FALSE, genes will be in default anti-alphabetical order)
	### rank.df: dataframe required if order.by.rank==TRUE, must provide bindev csv results here

	#error conditions
	if(order.by.rank==TRUE & is.null(rank.df)) {stop("order.by.rank = TRUE but no dataframe provided to pull rank information from")}
	if(length(grep("^ENSG",plot.genes))!=length(plot.genes)) {stop("ENSG names not detected in gene list provided. Concert from gene names to ensembl IDs")}
	
	#reminder messages of what mbp norm means for output
	if(norm.to.mbp==TRUE) {
		cat("norm.to.mbp=TRUE : avg expression values in raw counts divided by raw counts MBP\n")
		mbp.gene = rownames(spe)[rowData(spe)$gene_name=="MBP"]
		if(length(grep(mbp.gene, plot.genes))==0) plot.genes= c(plot.genes, mbp.gene)
	}

	if(norm.to.mbp==FALSE) cat("norm.to.mbp=FALSE : avg expression values in logNormCounts\n")

	#subset spe and build spe list
	spe2 = spe[plot.genes,]
	l4 = unique(spe2$sample_id)
	names(l4) = lapply(l4, function(x) unique(colData(spe2)[spe2$sample_id==x,"brain"]))
	l4 = lapply(l4, function(x) spe2[,colData(spe2)$sample_id==x])

	#build dotplot dataframe
	plot.genes.df = do.call(rbind, lapply(seq_along(l4), function(x) {
		if(norm.to.mbp==TRUE) {m1 = cbind("avg"=rowMeans(counts(l4[[x]])),
			"n"=rowSums(counts(l4[[x]])>0)/ncol(l4[[x]])) %>% as.data.frame() %>% tibble::rownames_to_column(var="gene")}
		if(norm.to.mbp==FALSE) {m1 = cbind("avg"=rowMeans(logcounts(l4[[x]])),
			"n"=rowSums(logcounts(l4[[x]])>0)/ncol(l4[[x]])) %>% as.data.frame() %>% tibble::rownames_to_column(var="gene")}
		
		m1$gene_name = rowData(l4[[x]])[m1$gene,"gene_name"]
		m1$brain = names(l4)[x]
		m1$position = unique(l4[[x]]$position)
		m1$slide = unique(l4[[x]]$slide)
		return(m1)
		})
	)
	
	#perform mbp norm if needed
	if(norm.to.mbp==TRUE) {
		mbp.df = filter(plot.genes.df, gene_name=="MBP")
		plot.genes.df = left_join(filter(plot.genes.df, gene_name!="MBP"), 
			mbp.df[,c("avg","brain","position","slide")], 
			by=c("brain", "position", "slide"), suffix=c("","_mbp")) %>% 
		mutate(expr_norm_mbp = avg/avg_mbp) %>% group_by(gene) %>% mutate(scaled_expr_norm_mbp= as.numeric(scale(expr_norm_mbp)))
	}

	#add in rank information if requestioned
	if(order.by.rank==TRUE | !is.null(rank.df)) {
		#create necessary rank dfs
		best.rank.slide.df = filter(rank.df, gene %in% plot.genes) %>% group_by(slide, gene, gene_name) %>%
			summarize(best.rank.slide=min(rank_brain), .groups="drop")
		best.rank.all.df = filter(rank.df, gene %in% plot.genes) %>% group_by(gene, gene_name) %>%
			summarize(best.rank.all = min(rank_brain), .groups="drop") %>%
			arrange(best.rank.all) %>% mutate(index=row_number(), i2=length(plot.genes)-index, ytext=paste(best.rank.all,gene_name, sep=" - ")) %>%
			arrange(i2) %>% mutate(ylabel=factor(i2, levels=i2,labels=ytext))
		#add rank information
		plot.genes.df <- left_join(plot.genes.df, best.rank.slide.df, by=c("slide","gene","gene_name")) %>%
			left_join(best.rank.all.df, by=c("gene","gene_name")) %>%
			mutate(xlabel=paste(position, brain))
	}
	return(plot.genes.df)
}
