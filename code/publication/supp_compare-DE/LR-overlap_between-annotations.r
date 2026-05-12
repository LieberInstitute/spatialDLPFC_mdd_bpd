setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

is_sig = function(x) x!="NS"

lr.sm = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv")
lr.se = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv")


compList = lapply(comparisons2, function(i) {
  cnames = grep(paste0(i,"_ttest"), colnames(lr.sm), value=T)
  clusnames = gsub(paste0("_",i,"_ttest"), "", cnames)
  
  tmp = cbind(select(lr.sm, gene_id, gene_name, "n_ttest_sig"=paste0("n_ttest_sig_",i)),
              lr.sm[,cnames]) %>%
    filter(n_ttest_sig>0) %>%
    mutate_at(cnames, is_sig)
  colnames(tmp)[4:ncol(tmp)] <- clusnames
  return(tmp)
})

compList2 = lapply(comparisons2, function(i) {
  cnames = grep(paste0(i,"_ttest"), colnames(lr.se), value=T)
  clusnames = gsub(paste0("_",i,"_ttest"), "", cnames)

  tmp = cbind(select(lr.se, gene_id, gene_name, "n_ttest_sig"=paste0("n_ttest_sig_",i)),
              lr.se[,cnames]) %>%
    filter(n_ttest_sig>0) %>%
    mutate_at(cnames, is_sig)
  colnames(tmp)[4:ncol(tmp)] <- clusnames
  return(tmp)
})

plot.df = do.call(rbind, lapply(comparisons2, function(i) {
	tmp_sm = compList[[i]]
	tmp_se = compList2[[i]]
	sm_names = colnames(tmp_sm)[4:ncol(tmp_sm)]
	se_names = colnames(tmp_se)[4:ncol(tmp_se)]
	out.df = do.call(rbind,	lapply(sm_names, function(j) {
                        t1 = tmp_sm[tmp_sm[[j]],"gene_name"]
			jcl = lapply(se_names, function(x) {
				q1 = tmp_se[tmp_se[[x]], "gene_name"]
				u1 = union(t1, q1)
				if(length(u1)==0) {
					out1 = data.frame("query"=x, jc=0, reference=j, n_union=0)
				} else {
					out1 = data.frame("query"=x, jc=length(intersect(t1, q1))/length(u1), reference=j, n_union=length(u1))
				}
				return(out1)
			})
			return(do.call(rbind, jcl))
                }))
	return(mutate(out.df, reference=factor(reference, levels=sm_names), query=factor(query, levels=rev(c(se_names, "self"))),
		comparison=i))
})) %>% mutate(comparison=factor(comparison, levels=comparisons2))

p1 <- ggplot(plot.df, aes(x=reference, y=query, fill=jc))+
  geom_tile(color="black", linewidth=.3)+
  scale_fill_gradientn(colors=c("white",RColorBrewer::brewer.pal(n=5,"Greys")[2:5],"black"))+
  facet_wrap(vars(comparison), ncol=6)+
  scale_x_discrete(labels=c("L1","L2","L3.4","L5","L6","WM"))+
  scale_y_discrete(labels=c("Olg","L6","L5","Inh","L4","L2.3","Ast","M.V"))+
  labs(x="PRECAST", y="Seurat", fill="Jaccard\nlrDEGs")+
  theme_minimal()+theme(#legend.position="none", 
			text=element_text(size=8), aspect.ratio=1,
                        panel.grid.major=element_blank(), panel.grid.minor=element_blank())

pdf(file="plots/publication/supp_compare-DE/LR-overlap_jaccard.pdf", width=6.5, height=2.5)
p1
dev.off()

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
