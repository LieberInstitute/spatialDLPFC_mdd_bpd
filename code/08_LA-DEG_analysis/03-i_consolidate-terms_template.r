library(dplyr)
library(igraph)
library(ggplot2)

set.seed(123)

source("code/08_LA-DEG_analysis/jc-igraph_functions.r")
source("code/08_LA-DEG_analysis/fgsea_functions.r")

#gmt_db = "Reactome"
#gmt_db = "GO-BP"
gmt_db = "GO-CC"

#res_file = "smoothed-k9-1663"
res_file = "seurat-pc30"

x="NTC.BPD_M"

source(paste0("code/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "_keylist.r"))

if(res_file=="seurat-pc30") keylist = keylist[["seurat"]]
if(res_file=="smoothed-k9-1663") keylist = keylist[["smoothed"]]

fgsea.results = readRDS(paste0("processed-data/08_LA-DEG_analysis/", res_file, "_", gmt_db, "_LA-fgsea-list.rda"))

ledge = tibble::deframe(fgsea.results[[x]][,c("pathway","leadingEdge")])

tmp2 = read.csv(paste0("processed-data/08_LA-DEG_analysis/", 
                       gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, 
                       "_neuronal-fgsea_filtered-LA-DEG-leadingeEdge.csv"))

ledge = ledge[unique(tmp2$term)]
jc.df = jaccardFromList(ledge)
e.df= filter(jc.df, jaccard>0)

keepterms = keylist[[gmt_db]]
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
#check for removed genes
filter(tmp2, gene_name %in% setdiff(tmp2$gene_name, tmp2_revised$gene_name))

#term to gene jaccard
outlist = term2GeneIGRAPH(tmp2_revised, layout_style="kk", text_title=paste(gsub("_"," ",x), gmt_db, res_file))

#check plot for more complicated networks
pdf(file=paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)), "_",
                gmt_db, "_", res_file, "_consolidated-terms.pdf"), width=12, height=12)
plot(outlist$igraph, layout=outlist$layout, vertex.label.family="sans", vertex.size=6, edges.curved=T,
     vertex.frame.width=0,
     main=outlist[["title_text"]])
legend("bottomleft", legend=c("Depleted","Enriched","Gene"),
       pch=16, pt.cex=1, cex=.7,
       col=c("#CFEBF7","#FFC0B5","grey85"))
dev.off()
