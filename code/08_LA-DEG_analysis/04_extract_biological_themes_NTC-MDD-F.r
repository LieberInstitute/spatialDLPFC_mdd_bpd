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

x="NTC.MDD_F"
target_group = unlist(strsplit(x, "_"))[[1]]
target_sex = unlist(strsplit(x, "_"))[[2]]

gmt_dbl = c("Reactome","GO-BP","GO-CC")

outlist1  = readRDS(paste0("processed-data/08_LA-DEG_analysis/", gsub("_","-", gsub("\\.","-",x)), "_seurat-pc30_consolidated-terms.rda"))
outlist2  = readRDS(paste0("processed-data/08_LA-DEG_analysis/", gsub("_","-", gsub("\\.","-",x)), "_smoothed-k9-1663_consolidated-terms.rda"))

keepterms = c(union(outlist1[[1]], outlist2[[1]]), union(outlist1[[2]], outlist2[[2]]), union(outlist1[[3]], outlist2[[3]]))
length(keepterms) #45

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

transmis = c("GABA","Synap")
termList[["syn.func"]] = setdiff(do.call(c, sapply(transmis, function(x) grep(x, keepterms, value=T))),
                                 "Synapse Assembly (GO:0007416)")
nrn.df = filter(all.terms, term %in% termList[["syn.func"]])
#length(unique(nrn.df$gene_name)) #10
geneList[["syn.func"]] = unique(nrn.df$gene_name)

neuro = c("Assembly","Dendrit","Projection","Neuronal","Calcium")
termList[["neuro"]] = do.call(c, sapply(neuro, function(x) grep(x, keepterms, value=T)))
nrn.df = filter(all.terms, term %in% termList[["neuro"]])
#length(unique(nrn.df$gene_name)) #13
geneList[["neuro"]] = unique(nrn.df$gene_name)

#termList[["neuro"]] = grep("Neuronal", keepterms, value=T, ignore.case = T)
#nrn.df = filter(all.terms, term %in% termList[["neuro"]])
##length(unique(nrn.df$gene_name)) #11
#geneList[["neuro"]] = unique(nrn.df$gene_name)

#calcium only had 3 genes and should correlate or be obvious if high enough fold change

mito = c("Mito")
termList[["mito"]] = do.call(c, lapply(mito, function(x) grep(x, keepterms, value=T)))
atp.df = filter(all.terms, term %in% termList[["mito"]])
#length(unique(atp.df$gene_name)) #9
geneList[["mito"]] = unique(atp.df$gene_name)

#glia = c("Drug","Collagen","Hemato")
#termList[["glia"]] = setdiff(do.call(c, lapply(glia, function(x) grep(x, keepterms, value=T))),
#                             "Collagen-Containing Extracellular Matrix (GO:0062023)")
#glia.df = filter(all.terms, term %in% termList[["glia"]])
##length(unique(glia.df$gene_name)) #17
#geneList[["glia"]] = unique(glia.df$gene_name)

#switch to enriched terms


vasc = c("Hemo","Platelet","Angio")
termList[["vasc"]] = do.call(c, sapply(vasc, function(x) grep(x, keepterms, value=T)))
vasc.df = filter(all.terms, term %in% termList[["vasc"]])
#table(vasc.df$NES) #all up, 32
#length(unique(vasc.df$gene_name)) #20
geneList[["vasc"]] = unique(vasc.df$gene_name)

#vasc2 = c("Endo","Epi")
#termList[["bbb"]] = do.call(c, lapply(vasc2, function(x) grep(x, keepterms, value=T)))
#vasc.df2 = filter(all.terms, term %in% termList[["bbb"]])
##table(vasc.df2$NES) #all up, 13
##length(unique(vasc.df2$gene_name)) #11
#geneList[["bbb"]] = unique(vasc.df2$gene_name)

#termList[["ecm"]] = "Collagen-Containing Extracellular Matrix (GO:0062023)"
#geneList[["ecm"]] = filter(all.terms, term=="Collagen-Containing Extracellular Matrix (GO:0062023)")$gene_name

#termList[["cytoskel"]] = "Cytoskeleton (GO:0005856)"
#geneList[["cytoskel"]] = filter(all.terms, term=="Cytoskeleton (GO:0005856)")$gene_name

#termList[["adhesion"]] = "Focal Adhesion (GO:0005925)"
#geneList[["adhesion"]] = filter(all.terms, term=="Focal Adhesion (GO:0005925)")$gene_name

inflamm = c("Cytokine","Inflamm","Neutrophil")
termList[["inflamm"]] = do.call(c, sapply(inflamm, function(x) grep(x, keepterms, value=T)))
inf.df = filter(all.terms, term %in% termList[["inflamm"]])
table(inf.df$NES) #all up, 94
length(unique(inf.df$gene_name)) #62
geneList[["inflamm"]] = unique(inf.df$gene_name)

#termList[["inflamm"]] = grep("Inflamm", keepterms, value=T)
#geneList[["inflamm"]] = unique(filter(all.terms, term %in% termList[["inflamm"]])$gene_name)

#termList[["neutrophil"]] = "Neutrophil Degranulation R-HSA-6798695"
#geneList[["neutrophil"]] = filter(all.terms, term=="Neutrophil Degranulation R-HSA-6798695")$gene_name

#termList[["cyto.prod"]] = "Regulation of Cytokine Production (GO:0001817)"
#geneList[["cyto.prod"]] = filter(all.terms, term=="Regulation of Cytokine Production (GO:0001817)")$gene_name

#termList[["cyto.resp"]] = setdiff(grep("Cytokine", keepterms, value=T), "Regulation of Cytokine Production (GO:0001817)")
#geneList[["cyto.resp"]] = unique(filter(all.terms, term %in% termList[["cyto.resp"]])$gene_name)

#termList[["complement"]] = "Complement Cascade R-HSA-166658"
#geneList[["complement"]] = filter(all.terms, term=="Complement Cascade R-HSA-166658")$gene_name

termList[["growth"]] = c("Negative Regulation of Growth (GO:0045926)",
                         "Negative Regulation of Cell Population Proliferation (GO:0008285)")
grow.df = filter(all.terms, term %in% termList[["growth"]])
geneList[["growth"]] = unique(grow.df$gene_name)

ribo = c("mRNA","P-body","NMD")
termList[["mRNA"]] = do.call(c, lapply(ribo, function(x) grep(x, keepterms, value=T)))
ribo.df = filter(all.terms, term %in% termList[["mRNA"]])
#table(ribo.df$NES) #all up, 50
#length(unique(ribo.df$gene_name)) #21
geneList[["mRNA"]] = unique(ribo.df$gene_name)

ribo = c("Ribo","Translation")
termList[["ribo"]] = do.call(c, lapply(ribo, function(x) grep(x, keepterms, value=T)))
ribo.df = filter(all.terms, term %in% termList[["ribo"]])
#table(ribo.df$NES) #all up, 50
#length(unique(ribo.df$gene_name)) #21
geneList[["ribo"]] = unique(ribo.df$gene_name)

termList[["starve"]] = c("Cellular Response To Starvation R-HSA-9711097","Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency R-HSA-9633012")
geneList[["starve"]] = unique(filter(all.terms, term %in% c("Cellular Response To Starvation R-HSA-9711097","Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency R-HSA-9633012"))$gene_name)
#Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency R…

termList[["apop"]] = grep("Apop", keepterms, value=T)
geneList[["apop"]] = unique(filter(all.terms, term %in% termList[["apop"]])$gene_name)

#filter(all.terms, !gene_name %in% unlist(geneList)) %>% arrange(gene_name)
#look up these genes to see if can place them in a theme

#remaining genes
cat("\n\nGenes that were not assigned to a theme:\n\n")
filter(all.terms[,1:4], !gene_name %in% unlist(geneList)) %>% arrange(gene_name)



#pull out top L-A DEGs with greatest logFC per theme
box_data = loadData()

top.genes = lapply(geneList, function(x) {
  g1 = filter(box_data[["LA_smoothed"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  g2 = filter(box_data[["LA_seurat"]], group==target_group, sex==target_sex, adj.P.Val<.05, gene_name %in% x) %>%
    slice_max(n=8, abs(logFC)) %>% pull(gene_name)
  intersect(g1, g2)
})

sapply(top.genes, length)


#plot themes
igraphList = lapply(names(geneList), function(y) {
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
spe_sm = box_data[["spe_smoothed"]][, box_data[["spe_smoothed"]]$condition %in% c("NTC","MDD") & box_data[["spe_smoothed"]]$sex=="F"]
spe_se = box_data[["spe_seurat"]][, box_data[["spe_seurat"]]$condition %in% c("NTC","MDD") & box_data[["spe_seurat"]]$sex=="F"]

setdiff(unique(all.terms$gene_name), rowData(spe_sm)$gene_name)
#S100A9, S100A8, and MT1A aren't present in smoothed pseudobulk so have to remove
setdiff(unique(all.terms$gene_name), rowData(spe_se)$gene_name)

hmp_data = formatData(setdiff(unique(all.terms$gene_name), c("MT1A","S100A9","S100A8")), spe_sm, spe_se)

#format for heatmap
d1 = tidyr::pivot_wider(hmp_data$boxplot.df, names_from="key_genes", values_from="logcounts")
#m1 = as.matrix(d1[,unique(all.terms$gene_name)])
m1 = as.matrix(d1[,6:ncol(d1)])
c1 = cor(m1, method="spearman")

d2 = distinct(ungroup(all.terms), gene_name, NES) %>%
  filter(!gene_name %in% c("MT1A","S100A9","S100A8"))

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
col.annot_d = col.annot[,c("NES","mito","neuro","syn.func")]
phmd = pheatmap(c1[dep.names,dep.names], annotation_col=col.annot_d[dep.names,], annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.MDD F L-A depleted genes\nTop logFC DEGs per theme annotated with black",
         silent=T)
#plot(phmd[[4]])

#plot enriched genes
en.names = rownames(col.annot)[col.annot$NES>0]
col.annot_e = col.annot[,c("NES", setdiff(names(geneList),c("mito","neuro","syn.func")))]
phme = pheatmap(c1[en.names,en.names], annotation_col=col.annot_e[en.names,], annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.MDD F L-A enriched genes\nTop logFC DEGs per theme annotated with black",
         silent=T)
#plot(phme[[4]])

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
spe_summ = spe_summ[,spe_summ$condition %in% c("NTC","MDD") & spe_summ$sex=="F"]
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
       y="GSEA term genes", title="NTC.MDD F L-A depleted terms (***top L-A DEGs per theme)")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        panel.background = element_rect(fill="#CFEBF7"),
                        panel.grid = element_line(color="#9FD7EF"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))
#p1

#enriched
pc30.df = dotplotDF(spe_summ, en.names, swap_rownames="gene_name", summarize_groups=T, 
                    cluster_labels="seurat_pc30", row_data=NULL) 
gene_order = rev(phme$tree_row$label[phme$tree_row$order])
pc30.df$gene_name_f = factor(pc30.df$gene_name, levels=gene_order,
                             labels=ifelse(gene_order %in% unlist(top.genes), paste0("***", gene_order), gene_order))
#dlpfc marker dotplot
#try splitting into two so that gene names are readable
p2 <- ggplot(filter(pc30.df, gene_name %in% phme$tree_row$label[phme$tree_row$order][1:grep("SCIN", phme$tree_row$label[phme$tree_row$order])]), 
	aes(x=factor(clusters, levels=seurat_levels), 
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.MDD F L-A enriched terms (***top L-A DEGs per theme)")+
  theme_minimal()+theme(axis.title.x=element_blank(), 
                        panel.background = element_rect(fill="#FFC0B5"),
                        panel.grid = element_line(color="#FFA190"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

p2.1 <- ggplot(filter(pc30.df, gene_name %in% phme$tree_row$label[phme$tree_row$order][(grep("SCIN", phme$tree_row$label[phme$tree_row$order])+1):length(phme$tree_row$label[phme$tree_row$order])]),
        aes(x=factor(clusters, levels=seurat_levels),
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.MDD F L-A enriched terms (***top L-A DEGs per theme)")+
  theme_minimal()+theme(axis.title.x=element_blank(),
                        panel.background = element_rect(fill="#FFC0B5"),
                        panel.grid = element_line(color="#FFA190"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))
#plot top genes expression with combo of violin and boxplots
## using custom function from 08_LA-DEG_analysis/plot-gex_functions.r
#depleted
plist1 = lapply(c("mito","neuro","syn.func"), function(x) {
  gl1 = top.genes[[x]]
  plist = lapply(gl1, function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
  plist = marrangeGrob(plist, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = paste(x, "top genes"))
  return(plist)
})
#enriched
plist2 = lapply(setdiff(names(top.genes), c("mito","neuro","syn.func")), function(x) {
  gl1 = top.genes[[x]]
  plist = lapply(gl1, function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
  plist = marrangeGrob(plist, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = paste(x, "top genes"))
  return(plist)
})

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
UpSetR::upset(UpSetR::fromList(geneList[c("mito","neuro","syn.func")]), nsets=3, text.scale=2, mb.ratio=c(.5,.5))
#UpSetR::upset(UpSetR::fromList(geneList[setdiff(names(geneList), c("energy","syn.struc","neuro","syn.func","glia"))]), 
#	nsets=length(geneList)-5, text.scale=2, mb.ratio=c(.5,.5))
plot(phmd[[4]])
p1
plist1
UpSetR::upset(UpSetR::fromList(geneList[setdiff(names(geneList), c("mito","neuro","syn.func"))]),
        nsets=length(geneList)-3, text.scale=2, mb.ratio=c(.5,.5))
plot(phme[[4]])
p2
p2.1
plist2
dev.off()
cat("\n\nPlots saved to:", paste0("plots/08_LA-DEG_analysis/", gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)),
                "_biological-themes.pdf"),"\n")

#pull strongly sig genes not in top.genes
### some genes were just barely not at logFC=.5 in only 1 so i lowered the threshold just a little but (genes like ADRA1D, CARNS1)
t1 = filter(box_data[["LA_smoothed"]], group==target_group, sex==target_sex, abs(logFC)>.3, adj.P.Val<.05)$gene_name
	#, !gene_name %in% unlist(top.genes))$gene_name
t2 = filter(box_data[["LA_seurat"]], group==target_group, sex==target_sex, abs(logFC)>.3, adj.P.Val<.05)$gene_name
	#, !gene_name %in% unlist(top.genes))$gene_name
new.genes = intersect(t1, t2)

#filter by prop n50
spe_save <- spe_summ
load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")
prop.m = t(assay(spe_summ, "logcounts.prop.detected")[rowData(spe_summ)$gene_name %in% new.genes,])
colnames(prop.m) = rowData(spe_summ)[colnames(prop.m), "gene_name"]

f2 <- function(x) as.numeric(median(x))
q2 = apply(prop.m, MARGIN=2, f2)
q2.pass = q2>.02
new.genes = names(q2.pass)[q2.pass]

spe_summ <- spe_save
hmp_data = formatData(union(unique(unlist(top.genes)), new.genes), spe_sm, spe_se)

#format for heatmap
d1 = tidyr::pivot_wider(hmp_data$boxplot.df, names_from="key_genes", values_from="logcounts")
m1 = as.matrix(d1[,union(unique(unlist(top.genes)), new.genes)])
c1 = cor(m1, method="spearman")

#d2 = distinct(ungroup(all.terms), gene_name, NES)
d2 = bind_rows(filter(box_data[["LA_seurat"]][,c("gene_name","logFC","group","sex")], group==target_group, sex==target_sex, gene_name %in% union(unlist(top.genes), new.genes)),
	filter(box_data[["LA_smoothed"]][,c("gene_name","logFC","group","sex")], group==target_group, sex==target_sex, gene_name %in% union(unlist(top.genes), new.genes))) %>%
	group_by(gene_name) %>% summarise("NES"=sign(mean(logFC)))

for(i in names(geneList)) {
  d2[[i]] = as.character(d2$gene_name %in% geneList[[i]])
}

#create and format heatmap columns annotations
col.annot = as.data.frame(d2[,-1])
rownames(col.annot) = d2$gene_name
col.annot$NES = as.character(col.annot$NES)

col.annot$logfc03 = as.character(rownames(col.annot) %in% new.genes)
col.annot$logfc03 = ifelse(rownames(col.annot) %in% unlist(top.genes) & col.annot$logfc03=="TRUE", "top", col.annot$logfc03)

for(i in names(top.genes)) {
  col.annot[top.genes[[i]],i] = "top"
}
#make logfc05 be bottom after NES
col.annot = col.annot[,c("NES","logfc03", names(top.genes))]

#select column annotation colors
annot_colors= list("NES"=c("-1"="skyblue","1"="tomato"),
	"logfc03"= c("FALSE"="white", "TRUE"="grey50", "top"="black"))
for(i in names(geneList)) {
  annot_colors[[i]] = c("FALSE"="white", "TRUE"="grey50", "top"="black")
}

#make c1 rownames have *** if top gene
rownames(c1) = ifelse(rownames(c1) %in% unlist(top.genes), paste0(rownames(c1), "***"), rownames(c1))

#plot heatmap
phm = pheatmap(c1, annotation_col=col.annot, annotation_colors = annot_colors,
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(9),
         breaks=seq(-1, 1, length.out=10), legend_breaks = seq(-1, 1, by=.5),
         angle_col=90, fontsize=7, treeheight_row = 12, treeheight_col = 12,
         clustering_method = "ward.D2", main="NTC.MDD F L-A top DEGs in term PLUS top DEGs (abs logFC>.3, zero prop n50= 2%) not in terms",
         silent=T)

#dotplot all genes with abs(logFC)>.5
gene_order = rev(intersect(phm$tree_col$label[phm$tree_col$order], new.genes))
pc30.df = dotplotDF(spe_summ, gene_order, swap_rownames="gene_name", summarize_groups=T,
                    cluster_labels="seurat_pc30", row_data=NULL)
#gene_order = rev(phm$tree_row$label[phm$tree_row$order])
pc30.df$gene_name_f = factor(pc30.df$gene_name, levels=gene_order,
                             labels=ifelse(gene_order %in% unlist(top.genes), paste0("***", gene_order), gene_order))

#add column with direction for background color
pc30.df2 = left_join(pc30.df, tibble::rownames_to_column(col.annot, var="gene_name"), by=c("gene_name")) %>%
	mutate(fill_color= factor(NES, levels=c("-1","1"), labels=c("#CFEBF7","#FFC0B5")))

#dlpfc marker dotplot
## split into two to make gene names legible
#sub1 = phm$tree_col$label[phm$tree_col$order][1:grep("DDIT4",phm$tree_col$label[phm$tree_col$order])]
sub1 = intersect(phm$tree_col$label[phm$tree_col$order], new.genes)[1:81]
#sub2 = phm$tree_col$label[phm$tree_col$order][(grep("DDIT4",phm$tree_col$label[phm$tree_col$order])+1):length(gene_order)]
sub2 = intersect(phm$tree_col$label[phm$tree_col$order], new.genes)[82:length(gene_order)]
p3 <- ggplot(filter(pc30.df, gene_name %in% sub1),
	aes(x=factor(clusters, levels=seurat_levels),
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=filter(pc30.df2, gene_name %in% sub1), 
	aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual(values=c("#CFEBF7"="#CFEBF7","#FFC0B5"="#FFC0B5"), guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.MDD F L-A DEGs with abs(logFC)>.3 and zero prop.>2% in half of capture areas\n(background indicates whether term was enriched/depleted)")+
  theme_minimal()+theme(axis.title.x=element_blank(),
                        #panel.background = element_rect(fill="palegoldenrod"),
                        #panel.grid = element_line(color="khaki3"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

p3.1 <- ggplot(filter(pc30.df, gene_name %in% sub2),
        aes(x=factor(clusters, levels=seurat_levels),
                          y=gene_name_f, color=mean_expr_scaled, size=prop_spots))+
  geom_tile(data=filter(pc30.df2, gene_name %in% sub2),
        aes(fill=fill_color, size=NA), alpha=.5, color="transparent", show.legend=c(size=FALSE))+
  scale_fill_manual(values=c("#CFEBF7"="#CFEBF7","#FFC0B5"="#FFC0B5"), guide="none")+
  geom_count()+scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=c("Micro\nVasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))+
  scale_size(range=c(1,6), limits=c(0,1), breaks=c(0,.5,1))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="GSEA term genes", title="NTC.MDD F L-A DEGs with abs(logFC)>.3 and zero prop.>2% in half of capture areas\n(background indicates whether term was enriched/depleted)")+
  theme_minimal()+theme(axis.title.x=element_blank(),
                        #panel.background = element_rect(fill="palegoldenrod"),
                        #panel.grid = element_line(color="khaki3"),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"),
                        axis.title.y=element_text(margin=margin(0,20,0,20,"pt")))

#gex plot all genes with abs(logFC)>.5
plist4 = lapply(intersect(phm$tree_col$label[phm$tree_col$order],new.genes), function(y) plotViolin(y, hmp_data$boxplot.df, box_data, color_by="cluster"))
plist4 = marrangeGrob(plist4, layout_matrix=matrix(c(1:8), ncol=2, nrow=4, byrow=T), top = "All genes with abs(logFC)>.3")

pdf(file=paste0("plots/08_LA-DEG_analysis/",gsub("\\.","-", gsub("_","-",x)),"/", gsub("\\.","-", gsub("_","-",x)), "_all-genes-logFC-03-zeroprop-02.pdf"), height=12, width=12)
plot(phm[[4]])
p3
p3.1
plist4
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
