library(dplyr)
library(igraph)
library(enrichR)

set.seed(123)

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

source("code/08_LR-DEG_analysis/fgsea_functions.r")
source("code/08_LR-DEG_analysis/jc-igraph_functions.r")
cpList = readRDS("plots/colorPalettes.rds")

gmt_db = "Reactome"
res_file = "smoothed-k9-1663"
clust_levels = c("L1","L2","L3.4","L5","L6","WM")
names(clust_levels) = clust_levels

#need full gmt list
if(gmt_db=="Reactome") gmt = .read_gmt("Reactome_2022")
if(gmt_db=="WikiPathways") gmt = .read_gmt("WikiPathways_2024_Human")

gmt.list = group_by(gmt, term) %>% summarise(gene=list(gene)) %>%
  tibble::deframe()


#load fgsea results
fgsea.results = readRDS(paste0("processed-data/08_LR-DEG_analysis/", res_file, "_", gmt_db, "_fgsea-list.rda"))

#load DEG lists
lrList = readRDS("processed-data/08_LR-DEG_analysis/LR-paired_rev-gene-input_logFC-0.3_lists.rds")
saveList <- readRDS("processed-data/07_dx_DE/LA-LR-overlap_rev-gene-input_lists.rds")

#in future will loop over all groups
groupList = names(fgsea.results)
names(groupList) <- groupList
#start with testing 1 group
x = "NTC.BPD_M"

#will loop overall clusters and make multi-page pdf
#start with testing 1 cluster
y = "L1"


#L1 JC igraph only
out1 = formatJaccardIGRAPH(fgsea.results, .group=x, .cluster=y)
outList = generateIGRAPH(out1)
#save plots
pdf(file="plots/08_LR-DEG_analysis/jc_test.pdf", height=8, width=8)
plot(simplify(outList[["igraph"]]), layout=outList[["layout"]], 
     edge.width=E(outList[["igraph"]])$jc*5, 
     vertex.label.family="sans", #vertex.label.font=2, 
     vertex.frame.color=NA, 
     main=paste(paste0(gsub("_", " ", x),":"), outList[["title_text"]]),
     sub=outList[["sub_text"]]
     )
legend("bottomleft", legend=c("Depleted","Depleted (with LR DEG)",
                                  "Enriched","Enriched (with LR DEG)"),
       pch=16, pt.cex=1, cex=.7,
       col=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
dev.off()


#also look at all terms from each cluster checked against one another
out2 = formatJaccardIGRAPH(fgsea.results, .group=x, .cluster=c("L1","L2","L3.4","L5","L6"))
outList2 = generateIGRAPH(out2)
#save plots
pdf(file="plots/08_LR-DEG_analysis/jc_test-all.pdf", height=12, width=12)
plot(simplify(outList2[["igraph"]]), layout=outList2[["layout"]], 
     edge.width=E(outList2[["igraph"]])$jc*5, 
     vertex.label.family="sans", #vertex.label.font=2, 
     vertex.frame.width=3,
     main=paste(paste0(gsub("_", " ", x),":"), outList2[["title_text"]]),
     sub=outList2[["sub_text"]]
)
legend("bottomleft", legend=c("Depleted","Depleted (with LR DEG)",
                              "Enriched","Enriched (with LR DEG)"),
       pch=16, pt.cex=1, cex=.7,
       col=c("#CFEBF7","skyblue","#FFC0B5","tomato"))
dev.off()

