setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})
set.seed(123)

prs.group = c("prsMDD","prsBPD","prsSCZ")

(results_set = "seurat-pc30-no-lowUMI")

if(results_set=="seurat-pc30-no-lowUMI") {
	results_set2 = "seurat-pc30"
} else {
	results_set2 = results_set
}

# load F test results
cat("\nLoad F-test results...\n")
cat("> PRS set\n")
prs.f = do.call(rbind, lapply(prs.group, function(x) {
  read.csv(paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set, "_rev-gene-input_F-test.csv")) %>%
    mutate(model=x)
})) %>% mutate(model=factor(model, levels=prs.group))
dim(prs.f)

cat("\n> dx-sex set\n")
la.f = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set2, "_rev-gene-input_F-test.csv"))
dim(la.f)

# load t test results
cat("\nLoad t-test results...\n")
cat("> PRS set\n")
prs.t = do.call(rbind, lapply(prs.group, function(x) {
  read.csv(paste0("processed-data/tmp_PRS_DE/layer-adjusted-pc3-age-nspots_", x, "_", results_set, "_rev-gene-input_moderated-t-test.csv"))
}))
dim(prs.t)

cat("\n> dx-sex	set\n")
la.t = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set2, "_rev-gene-input_moderated-t-test.csv"))
dim(la.t)

# scatter plot function
tScatterPlot <- function(GWAS_set, prs_t_df, dxsex_la_t_df, filter_F_sig=list(prs_f_df, dxsex_la_f_df),
                         plot.title=GWAS_set) {
  max.lim = round(max(max(abs(prs_t_df$t)),max(abs(dxsex_la_t_df$t))), 1)
  if(max.lim<max(max(abs(prs_t_df$t)),max(abs(dxsex_la_t_df$t)))) max.lim= max.lim+.1
  
  if(GWAS_set=="prsMDD") ntc.comp = "NTC.MDD"
  if(GWAS_set=="prsBPD") ntc.comp = "NTC.BPD"
  if(GWAS_set=="prsSCZ") ntc.comp = c("NTC.MDD","NTC.BPD")
  
  t.df = inner_join(filter(prs_t_df, group==GWAS_set, sex=="F"), 
                   filter(dxsex_la_t_df, coef %in% paste0("F_", c(ntc.comp, "MDD.BPD"))) %>% 
                     select(gene_name, gene_id, logFC, t, coef, sex, group),
                   by=c("gene_id", "gene_name","sex"), suffix=c("","_dx.sex")) %>%
    bind_rows(inner_join(filter(prs_t_df, group==GWAS_set, sex=="M"), 
                        filter(dxsex_la_t_df, coef %in% paste0("M_", c(ntc.comp, "MDD.BPD"))) %>% 
                          select(gene_name, gene_id, logFC, t, coef, sex, group),
                        by=c("gene_id", "gene_name", "sex"), suffix=c("","_dx.sex"))) %>%
    mutate(group_dx.sex=factor(group_dx.sex, levels=c(ntc.comp,"MDD.BPD")))
  
  if(!is.null(filter_F_sig)) {
    sig.f.genes = union(filter(filter_F_sig[[2]], adj.P.Val<.05)$gene_id,
                        filter(filter_F_sig[[1]], model==GWAS_set, adj.P.Val<.05)$gene_id)
    # to reflect truth that is subset to genes present in both pre-filtered pseudobulk sets
    sig.f.genes = intersect(sig.f.genes, t.df$gene_id)
    t.df = filter(t.df, gene_id %in% sig.f.genes)
    plot.title = paste0(plot.title, ": sig. F-test only (", length(sig.f.genes), " genes)")
  }
  
  r2.df = data.frame(group_dx.sex=rep(c(ntc.comp,"MDD.BPD"), 2),
                     sex=c(rep("F", length(ntc.comp)+1),rep("M", length(ntc.comp)+1)),
                     rsq=NA)
  for(i in 1:nrow(r2.df)) {
    check.r2 = filter(t.df, group_dx.sex==r2.df$group_dx.sex[[i]], sex==r2.df$sex[[i]])
    lm1 = lm(t ~ t_dx.sex, data=check.r2)
    r2.df[i,"rsq"] = signif(summary(lm1)$adj.r.squared, 3)
  }
  r2.df$group_dx.sex = factor(r2.df$group_dx.sex, levels=levels(t.df$group_dx.sex))
  
  ggplot(t.df, aes(x=t_dx.sex, y=t))+
    geom_point(size=.1)+geom_smooth(method="lm", se=F)+
    geom_text(data=r2.df, aes(x=-7, y=7, label=paste("R^2 =",rsq)), hjust=0)+
    ylim(-max.lim, max.lim)+xlim(-max.lim,max.lim)+
    facet_grid(cols=vars(group_dx.sex), rows=vars(sex))+
    labs(x="t statistic (dx-sex model)", y="t statistic (PRS*sex model)",
         title=plot.title)+
    theme_bw()+theme(aspect.ratio=1)
}

cat("\n\nPlot scatter...\n")
plist1 <- lapply(prs.group, function(x) tScatterPlot(GWAS_set=x, prs.t, la.t, filter_F_sig=NULL))

plist2 <- lapply(prs.group, function(x) tScatterPlot(GWAS_set=x, prs.t, la.t, filter_F_sig=list(prs.f, la.f)))

pdf(file=paste0("plots/tmp_PRS_DE/correlation-with-dxsex_", results_set, "_scatter.pdf"), height=7, width=6)
grid.arrange(rasterize(plist1[[1]], layer="Points", dpi=250), 
	rasterize(plist2[[1]], layer="Points", dpi=250), ncol=1)
grid.arrange(rasterize(plist1[[2]], layer="Points", dpi=250),
	rasterize(plist2[[2]], layer="Points", dpi=250), ncol=1)
grid.arrange(rasterize(plist1[[3]], layer="Points", dpi=250), 
	rasterize(plist2[[3]], layer="Points", dpi=250), ncol=1)
dev.off()

cat("\nScatter plots saved to:", paste0("plots/tmp_PRS_DE/correlation-with-dxsex_", results_set, "_scatter.pdf"),"\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

