setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(igraph)
	library(ggplot2)
})

set.seed(123)

source("code/08_LA-DEG_analysis/jc-igraph_functions.r")
source("code/08_LA-DEG_analysis/fgsea_functions.r")

x="NTC.MDD_F"

source(paste0("code/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_keylist.r"))

gmt_dbl = c("Reactome","GO-BP","GO-CC")

for(res_file in c("smoothed-k9-1663","seurat-pc30")) {

if(res_file=="seurat-pc30") klist = keylist[["seurat"]]
if(res_file=="smoothed-k9-1663") klist = keylist[["smoothed"]]

outlist <- list()
for(gmt_db in gmt_dbl) {
  cat(paste0("\n\n", gmt_db),"\n")
  
  fgsea.results = readRDS(paste0("processed-data/08_LA-DEG_analysis/", res_file, "_", gmt_db, "_LA-fgsea-list.rda"))
  
  ledge = tibble::deframe(fgsea.results[[x]][,c("pathway","leadingEdge")])
  
  tmp2 = read.csv(paste0("processed-data/08_LA-DEG_analysis/", 
                         gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, 
                         "_neuronal-fgsea_filtered-LA-DEG-leadingeEdge.csv"))
  
  ledge = ledge[unique(tmp2$term)]
  jc.df = jaccardFromList(ledge)
  e.df= filter(jc.df, jaccard>0)
  
  keepterms = klist[[gmt_db]]
  k_all = do.call(rbind, lapply(names(keepterms), function(x) {
    t1 = getTerms(x, e.df)
    return(data.frame("keyword"=rep(x, length(t1)), 
                      "term"=t1, "keep"= t1 %in% keepterms[[x]]))
  }))
  
  redundant.terms = filter(k_all, keep==F)$term
  keep.terms = filter(k_all, keep==T)$term
  #make sure not accidentally removing anything i want to keep
  redundant.terms = setdiff(redundant.terms,keep.terms)
  
  tmp2_revised = filter(tmp2, !term %in% redundant.terms)

  #extra for ntc mdd f
  if(x=="NTC.MDD_F" & gmt_db=="GO-BP") {
	extra.redundant = setdiff(tmp2_revised$term, keep.terms)
	tmp2_revised = filter(tmp2_revised, !term %in% extra.redundant)
  }
  #check for removed genes
  filter(tmp2, gene_name %in% setdiff(tmp2$gene_name, tmp2_revised$gene_name))
  
  #term to gene jaccard
  outlist[[gmt_db]] = unique(tmp2_revised$term)
}

saveRDS(outlist, paste0("processed-data/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_",
               res_file, "_consolidated-terms.rda"))
cat("\n\nSaved list of consolidated terms to:",paste0("processed-data/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_",
               res_file, "_consolidated-terms.rda"))
}

outlist1  = readRDS(paste0("processed-data/08_LA-DEG_analysis/",gsub("\\.","-", gsub("_","-",x)), "_seurat-pc30_consolidated-terms.rda"))
outlist2  = readRDS(paste0("processed-data/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_smoothed-k9-1663_consolidated-terms.rda"))

keepterms = c(union(outlist1[[1]], outlist2[[1]]), union(outlist1[[2]], outlist2[[2]]), union(outlist1[[3]], outlist2[[3]]))
length(keepterms)

# create plots
outlist = list()
for(gmt_db in gmt_dbl) {
  tmp2_sm = read.csv(paste0("processed-data/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), 
                            "_smoothed-k9-1663_", gmt_db, 
                            "_neuronal-fgsea_filtered-LA-DEG-leadingeEdge.csv"))
  tmp2_se = read.csv(paste0("processed-data/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), 
                            "_seurat-pc30_", gmt_db, 
                            "_neuronal-fgsea_filtered-LA-DEG-leadingeEdge.csv"))
  tmp3 = bind_rows(filter(tmp2_sm, term %in% keepterms) %>% mutate(NES=sign(NES)) %>% select(term, NES, source, gene_name, name),
                   filter(tmp2_se, term %in% keepterms) %>% mutate(NES=sign(NES)) %>% select(term, NES, source, gene_name, name))
  
  tmp3 = group_by(tmp3, term, NES, gene_name, name) %>% add_tally() %>%
    mutate(source=ifelse(n==2, "both", source)) %>% 
    distinct(term, NES, source, gene_name, name)
  
  outlist[[gmt_db]] = term2GeneIGRAPH(as.data.frame(tmp3), layout_style="kk", source=T, 
                             text_title=paste(gsub("_"," ",x), gmt_db, "both annotations"))
}
# save plots
pdf(file=paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)), 
                "_both-annotations_all-databases_consolidated-terms.pdf"), width=12, height=12)
for(i in outlist) {
  if(length(V(i$igraph))>120) {
    V(i$igraph)$size = ifelse(V(i$igraph)$type=="gene", 3, 6)
    V(i$igraph)$label.cex = .5
  } else {
    V(i$igraph)$size = 6
  }
  plot(i$igraph, layout=i$layout, vertex.label.family="sans", #vertex.size=6, 
       dges.curved=T,
       vertex.frame.width=0,
       edge.lty=as.character(factor(E(i$igraph)$source, 
                                    levels=c("both","seurat","smoothed"), 
                                    labels=c("solid","solid","dashed"))),
       edge.width=ifelse(E(i$igraph)$source=="both", 3, 1),
       main=i$title_text,
       sub="Thick, solid edges indicate sig. term-gene combo for both annotations. Thin, solid edges indicate sig. term-gene combo for Seurat label GSEA.\nThin, dashed edges indicate sig. term-gene combo for PRECAST (smoothed) GSEA.")
  legend("bottomleft", legend=c("Depleted","Enriched","Gene"),
         pch=16, pt.cex=1, cex=.7,
         col=c("#CFEBF7","#FFC0B5","grey85"))
}
dev.off()
cat("\n\nConsolidated term2Gene plots saved to:",paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)), 
	"_both-annotations_all-databases_consolidated-terms.pdf"))

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
