suppressPackageStartupMessages({
	library(dplyr)
	library(igraph)
	library(ggplot2)
	library(SpatialExperiment)
	library(pheatmap)
})

set.seed(123)

source("code/08_LA-DEG_analysis/jc-igraph_functions.r")
source("code/08_LA-DEG_analysis/plot-gex_functions.r")
# for dotplots
source("code/06_pseudobulk/custom_functions.r")


cpList <- readRDS("plots/colorPalettes.rds")

x="NTC.BPD_M"
target_group = unlist(strsplit(x, "_"))[[1]]
target_sex = unlist(strsplit(x, "_"))[[2]]

gmt_dbl = c("Reactome","GO-BP","GO-CC")

outlist1  = readRDS("processed-data/08_LA-DEG_analysis/NTC-BPD-M_seurat-pc30_consolidated-terms.rda")
outlist2  = readRDS("processed-data/08_LA-DEG_analysis/NTC-BPD-M_smoothed-k9-1663_consolidated-terms.rda")

keepterms = c(union(outlist1[[1]], outlist2[[1]]), union(outlist1[[2]], outlist2[[2]]), union(outlist1[[3]], outlist2[[3]]))
length(keepterms) #64

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
  outlist[[gmt_db]] = tmp3
  #outlist[[gmt_db]] = term2GeneIGRAPH(as.data.frame(tmp3), layout_style="kk", source=T, 
  #                           text_title=paste(gsub("_"," ",x), gmt_db, "both annotations"))
}

all.terms = do.call(rbind, outlist)

termList = list()
geneList = list()

cellular.upkeep = c("Translation","Electron","TCA","ATP","Respirat","Mitochondrial")
termList[["upkeep"]] = do.call(c, sapply(cellular.upkeep, function(x) grep(x, keepterms, value=T)))
cup.df = filter(all.terms, term %in% termList[["upkeep"]])
#table(cup.df$NES) #all down, 126
#length(unique(cup.df$gene_name)) #44
geneList[["upkeep"]] = unique(cup.df$gene_name)

proinfl = c("Cytokine","JNK","Inflammat","Apop")
termList[["proinfl"]] = setdiff(do.call(c, sapply(proinfl, function(x) grep(x, keepterms, value=T))),
                  c("Cell Recruitment (Pro-Inflammatory Response) R-HSA-9664424",
                    "Regulation of Apoptotic Signaling Pathway (GO:2001233)"))
pro.df = filter(all.terms, term %in% termList[["proinfl"]])
#table(pro.df$NES) #all up, 39
#length(unique(pro.df$gene_name)) #28
geneList[["proinfl"]] = unique(pro.df$gene_name)

nuclear = c("RNA","Chromatin","Nucle")
termList[["nuclear"]] = do.call(c, sapply(nuclear, function(x) grep(x, keepterms, value=T)))
nuc.df = filter(all.terms, term %in% termList[["nuclear"]])
#table(nuc.df$NES) #all up, 24
#length(unique(nuc.df$gene_name)) #21
geneList[["nuclear"]] = unique(nuc.df$gene_name)

vascular = c("Platelet","Endothelial","Epithelial","Blood","Mechanical","Artery")
termList[["vascular"]] = do.call(c, sapply(vascular, function(x) grep(x, keepterms, value=T)))
vasc.df = filter(all.terms, term %in% termList[["vascular"]] )
#table(vasc.df$NES) #5 down, 14 up
#length(unique(vasc.df$gene_name)) #17
geneList[["vascular"]] = unique(vasc.df$gene_name)


immune = c("Complement","MHC","Microglia","Phago","Immune","Chemo","Pro-Infl","Secretory")
termList[["immune"]]  = setdiff(do.call(c, sapply(immune, function(x) grep(x, keepterms, value=T))),
                 c("Cytokine Signaling In Immune System R-HSA-1280215",
                   "Secretory Granule Lumen (GO:0034774)"))
imm.df = filter(all.terms, term %in% termList[["immune"]])
#table(imm.df$NES) #all down, 45
#length(unique(imm.df$gene_name)) #18
geneList[["immune"]] = unique(imm.df$gene_name)

#temporary one for GPCR
gpcr = c("GPCR","G Alpha")
termList[["gpcr"]] = do.call(c, sapply(gpcr, function(x) grep(x, keepterms, value=T)))
gpcr.df = filter(all.terms, term %in% termList[["gpcr"]])
#table(gpcr.df$NES) #all down, 33
geneList[["gpcr"]] = unique(gpcr.df$gene_name)
#GPCR to immune 
#geneList[["immune"]] = c(geneList[["immune"]],
#                         c("CX3CR1","P2RY13","P2RY12","C3","C3AR1","RGS10","RGS5","PRKCB","LPAR6"))

#GPCR to GABA
geneList[["gaba"]] = c("SST","CORT","PNOC","CRH","TAC1")
gaba = grep("GABA", keepterms, value=T)
gaba.df = filter(all.terms, term %in% gaba)
termList[["gaba"]] = c("GPCR Ligand Binding R-HSA-500792", gaba)
geneList[["gaba"]] = c(geneList[["gaba"]], unique(gaba.df$gene_name))


#filter(all.terms, !gene_name %in% unlist(geneList)) %>% arrange(gene_name)
#look up these genes to see if can place them in a theme

geneList[["vascular"]] = c(geneList[["vascular"]], "ADM", "FOLR2")
geneList[["nuclear"]] = c(geneList[["nuclear"]], "BTG1")
geneList[["upkeep"]] = c(geneList[["upkeep"]], "CHCHD2", "TOMM7")
geneList[["gaba"]] = c(geneList[["gaba"]], "RELN")

#remaining genes
cat("\n\nGenes that were not assigned to a theme:\n\n")
filter(all.terms[,1:4], !gene_name %in% unlist(geneList)) %>% arrange(gene_name)
#other neuronal: ARF1, CP, GRIK1
#not included in any theme: CRISPLD2, PYGB, SMAD1, TGFBR2, UFC1, ZMIZ1, ZNF395


#pull out top L-A DEGs with greatest logFC per theme
box_data = loadData()

top.genes = lapply(geneList, function(x) {
  g1 = filter(box_data[["LA_smoothed"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  g2 = filter(box_data[["LA_seurat"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  intersect(g1, g2)
})

#sapply(top.genes, length)


#plot themes
igraphList = lapply(c("upkeep","gaba","gpcr","immune","vascular","proinfl","nuclear"), function(y) {
  t1 = filter(all.terms, gene_name %in% geneList[[y]])
  i1 = term2GeneIGRAPH(as.data.frame(t1), layout_style="kk",
                  text_title=paste(gsub("_"," ",x), y, "theme"))
  #make top DEG darker
  V(i1$igraph)$color = ifelse(V(i1$igraph)$name %in% top.genes[[y]], "grey50", V(i1$igraph)$color)
  #make theme terms darker
  mod.color = V(i1$igraph)$name %in% termList[[y]]
  V(i1$igraph)$color[mod.color] = ifelse(V(i1$igraph)$color[mod.color]=="#CFEBF7", "skyblue", "tomato")
  
  return(i1)
})



# plot gene-gene pseudobulk sample level correlations in heatmap form
spe_sm = box_data[["spe_smoothed"]][, box_data[["spe_smoothed"]]$condition %in% c("NTC","BPD") & box_data[["spe_smoothed"]]$sex=="M"]
spe_se = box_data[["spe_seurat"]][, box_data[["spe_seurat"]]$condition %in% c("NTC","BPD") & box_data[["spe_seurat"]]$sex=="M"]
hmp_data = formatData(unique(all.terms$gene_name), spe_sm, spe_se)

#format for heatmap
d1 = tidyr::pivot_wider(hmp_data$boxplot.df, names_from="key_genes", values_from="logcounts")
m1 = as.matrix(d1[,unique(all.terms$gene_name)])
c1 = cor(m1, method="spearman")

d2 = distinct(ungroup(all.terms), gene_name, NES)
for(i in names(geneList)) {
  d2[[i]] = as.character(d2$gene_name %in% geneList[[i]])
}

#create and format heatmap columns annotations
col.annot = as.data.frame(d2[,-1])
rownames(col.annot) = d2$gene_name
col.annot$NES = as.character(col.annot$NES)

for(i in names(top.genes)) {
  col.annot[top.genes[[i]],i] = "top"
}

#select column annotation colors
annot_colors= list("NES"=c("-1"="skyblue","1"="tomato"))
for(i in names(geneList)) {
  annot_colors[[i]] = c("FALSE"="white", "TRUE"="grey50", "top"="black")
}

#plot depleted genes
dep.names = rownames(col.annot)[col.annot$NES<0]
phmd = pheatmap(c1[dep.names,dep.names], annotation_col=col.annot[dep.names,], annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.BPD M L-A depleted genes\nTop logFC DEGs per theme annotated with black",
         silent=T)

#plot enriched genes
en.names = rownames(col.annot)[col.annot$NES>0]
phme = pheatmap(c1[en.names,en.names], annotation_col=col.annot[en.names,], annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.BPD M L-A enriched genes\nTop logFC DEGs per theme annotated with black",
         silent=T)

# plot dotplots to look at layer distribution of theme genes
## using custom function from 06_pseudobulk/custom_functions.r

#load in sce for heatmap/dotplots (decided just to use seurat pc30)
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

cat("\nSubset spe object for dotplot to only NTC M and BPD M samples:\n")
spe_summ = spe_summ[,spe_summ$condition %in% c("NTC","BPD") & spe_summ$sex=="M"]
dim(spe_summ)
#depleted
pc30.df = dotplotDF(spe_summ, dep.names, swap_rownames="gene_name", summarize_groups=T, 
                    cluster_labels="seurat_pc30", row_data=NULL) 
gene_order = rev(phmd$tree_row$label[phmd$tree_row$order])
pc30.df$gene_name_f = factor(pc30.df$gene_name, levels=gene_order,
                             labels=ifelse(gene_order %in% unlist(top.genes), paste0("***", gene_order), gene_order))

#dlpfc marker dotplot
p1 <- ggplot(pc30.df, aes(x=factor(clusters, levels=seurat_levels), 
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.BPD M L-A depleted terms (***top L-A DEGs per theme)")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        panel.background = element_rect(fill="#CFEBF7"),
                        panel.grid = element_line(color="#9FD7EF"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

#enriched
pc30.df = dotplotDF(spe_summ, en.names, swap_rownames="gene_name", summarize_groups=T, 
                    cluster_labels="seurat_pc30", row_data=NULL) 
gene_order = rev(phme$tree_row$label[phme$tree_row$order])
pc30.df$gene_name_f = factor(pc30.df$gene_name, levels=gene_order,
                             labels=ifelse(gene_order %in% unlist(top.genes), paste0("***", gene_order), gene_order))
#dlpfc marker dotplot
p2 <- ggplot(pc30.df, aes(x=factor(clusters, levels=seurat_levels), 
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.BPD M L-A enriched terms (***top L-A DEGs per theme)")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        panel.background = element_rect(fill="#FFC0B5"),
                        panel.grid = element_line(color="#FFA190"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

#plot top genes expression with combo of violin and boxplots
## using custom function from 08_LA-DEG_analysis/plot-gex_functions.r
#depleted
plist1 = lapply(c("upkeep","gaba","gpcr","immune"), function(x) {
  gl1 = top.genes[[x]]
  plist = lapply(gl1, function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
  plist = marrangeGrob(plist, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = paste(x, "top genes"))
  return(plist)
})
#enriched
plist2 = lapply(c("proinfl","nuclear","vascular"), function(x) {
  gl1 = top.genes[[x]]
  plist = lapply(gl1, function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
  plist = marrangeGrob(plist, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = paste(x, "top genes"))
  return(plist)
})

#now do a revised version of the enriched guys, splitting by glial and neuronal
gene_order = phme$tree_row$label[phme$tree_row$order]
en.names1 = gene_order[1:grep("CRISPLD2", gene_order)]
en.names2 = gene_order[(grep("CRISPLD2", gene_order)+1):length(gene_order)]

top.genes2 = lapply(list(en.names1, en.names2), function(x) {
  g1 = filter(box_data[["LA_smoothed"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  g2 = filter(box_data[["LA_seurat"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  intersect(g1, g2)
})

t1 = filter(all.terms, gene_name %in% en.names1)
i1 = term2GeneIGRAPH(as.data.frame(t1), layout_style="kk",
                     text_title=paste(gsub("_"," ",x), "enriched (glial)"))
#make top DEG darker
V(i1$igraph)$color = ifelse(V(i1$igraph)$name %in% top.genes2[[1]], "grey50", V(i1$igraph)$color)
#make theme terms darker
mod.color = V(i1$igraph)$name %in% termList[["proinfl"]]
V(i1$igraph)$color[mod.color] = ifelse(V(i1$igraph)$color[mod.color]=="#CFEBF7", "skyblue", "tomato")


t2 = filter(all.terms, gene_name %in% en.names2)
i2 = term2GeneIGRAPH(as.data.frame(t2), layout_style="kk",
                     text_title=paste(gsub("_"," ",x), "enriched (neuronal)"))
#make top DEG darker
V(i2$igraph)$color = ifelse(V(i2$igraph)$name %in% top.genes2[[2]], "grey50", V(i2$igraph)$color)
#make theme terms darker
mod.color = V(i2$igraph)$name %in% termList[["proinfl"]]
V(i2$igraph)$color[mod.color] = ifelse(V(i2$igraph)$color[mod.color]=="#CFEBF7", "skyblue", "tomato")

#violin for revised neuronal
plist3 = lapply(top.genes2[[1]], function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
plist3 = marrangeGrob(plist3, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = "Glial enriched top genes")

plist4 = lapply(top.genes2[[2]], function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
plist4 = marrangeGrob(plist4, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = "Neuronal enriched top genes")

# save plots
pdf(file=paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)), 
                "_biological-themes.pdf"), width=12, height=12)
for(i in igraphList) {
  if(length(V(i$igraph))>120) {
    V(i$igraph)$size = ifelse(V(i$igraph)$type=="gene", 3, 6)
    V(i$igraph)$label.cex = .5
  } else {
    V(i$igraph)$size = 6
  }
  plot(i$igraph, layout=i$layout, vertex.label.family="sans", #vertex.size=6, 
       vertex.frame.width=0,
       main=i$title_text,
       sub="Reactome, GO-BP, and GO-CC results combined to produce theme")
  legend("bottomleft", legend=c("Depleted term","Depleted theme term","Enriched term","Enriched theme term","Gene","Top L-A DEG"),
         pch=16, pt.cex=1, cex=.7,
         col=c("#CFEBF7","skyblue","#FFC0B5","tomato","grey85","grey50"))
}
UpSetR::upset(UpSetR::fromList(geneList), nsets=length(geneList), text.scale=2, mb.ratio=c(.5,.5))
plot(phmd[[4]])
p1
plist1
plot(phme[[4]])
p2
plist2
#now plot revised ones
for(i in list(i1, i2)) {
  if(length(V(i$igraph))>120) {
    V(i$igraph)$size = ifelse(V(i$igraph)$type=="gene", 3, 6)
    V(i$igraph)$label.cex = .5
  } else {
    V(i$igraph)$size = 6
  }
  plot(i$igraph, layout=i$layout, vertex.label.family="sans", #vertex.size=6, 
       vertex.frame.width=0,
       main=i$title_text,
       sub="Reactome, GO-BP, and GO-CC results combined to produce theme")
  legend("bottomleft", legend=c("Depleted term","Enriched term","Enriched proinfl term","Gene","Top L-A DEG"),
         pch=16, pt.cex=1, cex=.7,
         col=c("#CFEBF7","skyblue","#FFC0B5","tomato","grey85","grey50"))
}
plist3
plist4
dev.off()
cat("\n\nPlots saved to:", paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)),
                "_biological-themes.pdf"),"\n")

#test additional plot
#pull strongly sig genes not in top.genes
t1 = filter(box_data[["LA_smoothed"]], group=="NTC.BPD", sex=="M", abs(logFC)>.5, adj.P.Val<.05)$gene_name
	#, !gene_name %in% unlist(top.genes))$gene_name
t2 = filter(box_data[["LA_seurat"]], group=="NTC.BPD", sex=="M", abs(logFC)>.5, adj.P.Val<.05)$gene_name
	#, !gene_name %in% unlist(top.genes))$gene_name
new.genes = intersect(t1, t2)
hmp_data = formatData(union(unique(unlist(top.genes)), new.genes), spe_sm, spe_se)

#format for heatmap
d1 = tidyr::pivot_wider(hmp_data$boxplot.df, names_from="key_genes", values_from="logcounts")
m1 = as.matrix(d1[,union(unique(unlist(top.genes)), new.genes)])
c1 = cor(m1, method="spearman")

#d2 = distinct(ungroup(all.terms), gene_name, NES)
d2 = bind_rows(filter(box_data[["LA_seurat"]][,c("gene_name","logFC","group","sex")], group=="NTC.BPD", sex=="M", gene_name %in% union(unlist(top.genes), new.genes)),
	filter(box_data[["LA_smoothed"]][,c("gene_name","logFC","group","sex")], group=="NTC.BPD", sex=="M", gene_name %in% union(unlist(top.genes), new.genes))) %>%
	group_by(gene_name) %>% summarise("NES"=sign(mean(logFC)))

for(i in names(geneList)) {
  d2[[i]] = as.character(d2$gene_name %in% geneList[[i]])
}

#create and format heatmap columns annotations
col.annot = as.data.frame(d2[,-1])
rownames(col.annot) = d2$gene_name
col.annot$NES = as.character(col.annot$NES)

col.annot$logfc05 = as.character(rownames(col.annot) %in% new.genes)
col.annot$logfc05 = ifelse(rownames(col.annot) %in% unlist(top.genes) & col.annot$logfc05=="TRUE", "top", col.annot$logfc05)

for(i in names(top.genes)) {
  col.annot[top.genes[[i]],i] = "top"
}
#make logfc05 be bottom after NES
col.annot = col.annot[,c("NES","logfc05", names(top.genes))]

#select column annotation colors
annot_colors= list("NES"=c("-1"="skyblue","1"="tomato"),
	"logfc05"= c("FALSE"="white", "TRUE"="grey50", "top"="black"))
for(i in names(geneList)) {
  annot_colors[[i]] = c("FALSE"="white", "TRUE"="grey50", "top"="black")
}
#annot_colors[["new.genes"]] = c("FALSE"="white", "TRUE"="grey50")

#make c1 rownames have *** if top gene
rownames(c1) = ifelse(rownames(c1) %in% unlist(top.genes), paste0(rownames(c1), "***"), rownames(c1))

#plot heatmap
phm = pheatmap(c1, annotation_col=col.annot, annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.BPD M L-A top DEGs in term supplemented with top DEGs (abs logFC>.5) not in terms",
         silent=T)

#dotplot all genes with abs(logFC)>.5
gene_order = rev(phm$tree_col$label[phm$tree_col$order])
pc30.df = dotplotDF(spe_summ, gene_order, swap_rownames="gene_name", summarize_groups=T,
                    cluster_labels="seurat_pc30", row_data=NULL)
#gene_order = rev(phm$tree_row$label[phm$tree_row$order])
pc30.df$gene_name_f = factor(pc30.df$gene_name, levels=gene_order,
                             labels=ifelse(gene_order %in% unlist(top.genes), paste0("***", gene_order), gene_order))

#add column with direction for background color
pc30.df2 = left_join(pc30.df, tibble::rownames_to_column(col.annot, var="gene_name"), by=c("gene_name")) %>%
	mutate(fill_color= factor(NES, levels=c("-1","1"), labels=c("#CFEBF7","#FFC0B5")))

#dlpfc marker dotplot
p3 <- ggplot(pc30.df, aes(x=factor(clusters, levels=seurat_levels),
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=pc30.df2, aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual(values=c("#CFEBF7","#FFC0B5"), guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.BPD M L-A DEGs with abs(logFC)>.5 (background indicates whether term was enriched/depleted)")+
  theme_minimal()+theme(axis.title.x=element_blank(),
                        #panel.background = element_rect(fill="palegoldenrod"),
                        #panel.grid = element_line(color="khaki3"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

#gex plot all genes with abs(logFC)>.5
plist4 = lapply(phm$tree_col$label[phm$tree_col$order], function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
plist4 = marrangeGrob(plist4, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = "All genes with abs(logFC)>.5")

pdf(file="plots/08_LA-DEG_analysis/NTC-BPD-M/NTC-BPD-M_all-genes-logFC-05.pdf", height=12, width=12)
plot(phm[[4]])
p3
plist4
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
