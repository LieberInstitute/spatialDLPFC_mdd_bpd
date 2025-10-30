library(dplyr)
library(igraph)
library(pheatmap)
library(ggplot2)

set.seed(123)
                
source("code/08_LR-DEG_analysis/jc-igraph_functions.r")

cpList = readRDS("plots/colorPalettes.rds")
cpList$seurat.bright = cpList$transfer.bright
cpList$seurat.light = cpList$transfer.light

gmt_db = "Reactome"
#gmt_db = "GO-BP"
#gmt_db = "GO-CC"

#res_file = "smoothed-k9-1663"

res_file = "seurat-pc30"

x="NTC.BPD_M"

tmp2 = read.csv(paste0("processed-data/08_LR-DEG_analysis/", 
                       gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, 
                       "_neuronal-fgsea_filtered-DEG-leadingeEdge.csv"))
out2 = readRDS(paste0("processed-data/08_LR-DEG_analysis/", 
                      gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, 
                      "_neuronal-fgsea_jc-igraph.rda"))


keywords = c("GPCR","Citric Acid","Translation","Complement","Interleukin")
names(keywords) <- keywords


getTerms <- function(keyword, edges_DF) {
  keyword_terms = union(unique(grep(keyword, edges_DF$reference, ignore.case=T,value=T)),
                        unique(grep(keyword, edges_DF$query, ignore.case=T,value=T)))
  
  keyword_subset = filter(edges_DF, reference %in% keyword_terms | query %in% keyword_terms)
  #plot(density(keyword_subset[,"jaccard"]))
  keyword_subset2 = filter(keyword_subset, jaccard>.2)
  unique.terms = unique(c(keyword_subset2$reference, keyword_subset2$query, keyword_terms))
  unique.terms = unique(unlist(lapply(strsplit(unique.terms, " "), function(x) paste(x[-1], collapse=" "))))
  return(unique.terms)
}

keyterms = lapply(keywords, getTerms, edges_DF=out2$edges)


i1 = lapply(keyterms, term2GeneIGRAPH, leadingEdge_DF=tmp2)

pdf(file=paste0("plots/08_LR-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)), "/",
                gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_keyword-igraph.pdf"), height=8, width=8)
for(i in names(i1)) {
  grid.size =factor(length(unique(i1[[i]]$edges$cluster)), levels=1:6, 
                    labels=c("single","small","medium","medium","large","large"))
  if(grid.size=="single") par(mfrow=c(1,1))
  if(grid.size=="small") par(mfrow=c(1,2))
  if(grid.size=="medium") par(mfrow=c(2,2)) 
  if(grid.size=="large") par(mfrow=c(2,3))
  for(j in unique(i1[[i]]$edges$cluster)) {
    plot(subgraph_from_edges(i1[[i]]$igraph, which(E(i1[[i]]$igraph)$cluster==j), delete.vertices = F), 
         layout=i1[[i]]$layout, 
         vertex.label.family="sans", size=8, edges.curved=F, 
         vertex.frame.width=0, asp=0, ylim=c(-1.1,1.1),
         main=paste(gsub("_"," ",x), paste0(gmt_db, ":"), i, paste0("(",j,")")))
  }
}
dev.off()

#decided terms
decidedterms <- list(c("GPCR Ligand Binding R-HSA-500792","GPCR Downstream Signaling R-HSA-388396"),
                     c("Respiratory Electron Transport R-HSA-611105"),
                     c("Translation R-HSA-72766"),
                     c("Complement Cascade R-HSA-166658"),
                     c("Cytokine Signaling In Immune System R-HSA-1280215"))
names(decidedterms) <- names(keyterms)

k_all = do.call(rbind, lapply(keywords, function(x) data.frame("keyword"=rep(x, length(keyterms[[x]])), 
                                                               "term"=keyterms[[x]], "keep"= keyterms[[x]] %in% decidedterms[[x]]))
)


redundant.terms = filter(k_all, keep==F)$term
length(redundant.terms) #17
keep.terms = filter(k_all, keep==T)$term
#make sure not accidentally removing anything i want to keep
redundant.terms = setdiff(redundant.terms,keep.terms)
length(redundant.terms) #17

tmp2_revised = filter(tmp2, !term %in% redundant.terms)
filter(tmp2, gene_name %in% setdiff(tmp2$gene_name, tmp2_revised$gene_name))
#JAK3

extra.terms = c("Platelet Degranulation R-HSA-114608")

redundant.terms = c(redundant.terms, extra.terms)
length(redundant.terms) #18
tmp2_revised = filter(tmp2, !term %in% redundant.terms)
filter(tmp2, gene_name %in% setdiff(tmp2$gene_name, tmp2_revised$gene_name))
#JAK3



sample(1:10,4)
dim(out2$nodes) #70
n.df = out2$nodes
n.df$term2 = unlist(lapply(strsplit(n.df$term, " "), function(x) paste(x[-1], collapse=" ")))
n.df = filter(n.df, !term2 %in% redundant.terms)
dim(n.df) #31 rows

dim(out2$edges) #659 rows
e.df = out2$edges
e.df$ref2 = unlist(lapply(strsplit(e.df$reference, " "), function(x) paste(x[-1], collapse=" ")))
e.df = filter(e.df, !ref2 %in% redundant.terms)
e.df$query2 = unlist(lapply(strsplit(e.df$query, " "), function(x) paste(x[-1], collapse=" ")))
e.df = filter(e.df, !query2 %in% redundant.terms)
dim(e.df) #102

out2_revised = list("nodes"=n.df[,1:8], "edges"=e.df[,1:5])


outList2 = generateIGRAPH(out2_revised)

i2 = term2GeneIGRAPH(unique(tmp2_revised$term), tmp2_revised)

pdf(file=paste0("plots/08_LR-DEG_analysis/",gsub("\\.","-", gsub("_","-",x)), "/",
                gsub("\\.","-", gsub("_","-",x)), "_", res_file, "_", gmt_db, "_rev-igraph.pdf"), height=12, width=12)
plot(simplify(outList2[["igraph"]]), layout=outList2[["layout"]], 
     edge.width=E(outList2[["igraph"]])$jc*5, 
     vertex.label.family="sans", #vertex.label.font=2, 
     vertex.frame.width=3,
     main=paste(paste0(gsub("_", " ", x)," ", gmt_db, ":"), outList2[["title_text"]]),
     sub=outList2[["sub_text"]]
)
legend("bottomleft", legend=c("Depleted", "Enriched"),
       pch=16, pt.cex=1, cex=.7,
       col=c("skyblue","tomato"))
set.seed(123)
plot(i2$igraph,
     layout=i2$layout, 
     #layout=layout.fruchterman.reingold(i2$igraph), #can use instead of fr layout if too crowded
     vertex.label.family="sans", vertex.size=6, edges.curved=T, 
     vertex.frame.width=0, #asp=0, ylim=c(-1.1,1.1),
     main=paste(gsub("_"," ",x), paste0(gmt_db, ":"), factor(unlist(strsplit(res_file, "-"))[[1]], levels=c("smoothed","seurat"),
                                                     labels=c("PRECAST (smoothed)","Seurat labels"))),
     sub="Thin, light edges indicate L-A DEG. Darker, thicker edges indicate L-R DEG.\nDark grey nodes (if present) are only L-R DEGs and are not L-A significant.")
legend("bottomleft", 
       legend=names(cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]]),
       lty=1, lwd=2, cex=.7,
       col=cpList[[paste0(unlist(strsplit(res_file, "-"))[[1]],".bright")]])
dev.off()
