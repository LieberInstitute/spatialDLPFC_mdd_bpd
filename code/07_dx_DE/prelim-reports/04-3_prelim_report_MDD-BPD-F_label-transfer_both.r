setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scater)
	library(dplyr)
	library(ggplot2)
	library(clusterProfiler)
	library(enrichR)
	library(enrichplot)
	library(pheatmap)
	library(gridExtra)
	library(ggrastr)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
cat("\nLoading label transfer results...\n")

# load layer-adjusted results
adj.results = read.csv("processed-data/07_dx_DE/layer-adjusted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased", "increased")))

sig.la = filter(adj.results, sex=="F", group=="MDD.BPD", adj.P.Val<.01)


# load layer-restricted results
restr.results <- read.csv("processed-data/07_dx_DE/layer-restricted-age_seurat-pc30-no-lowUMI_compiled-results.csv", row.names=1) %>%
  mutate(sex= factor(sex, levels=c("F","M")),
         group= factor(group, levels=c("NTC.MDD","NTC.BPD","MDD.BPD")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb"))
         )

sig.lr = bind_rows(filter(restr.results, cluster!="Oligo", adj.P.Val<.01),
                   filter(restr.results, cluster=="Oligo", sex=="M", adj.P.Val<.01),
                   filter(restr.results, cluster=="Oligo", sex=="F", adj.P.Val<.0001)) %>%
  mutate(dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased","increased"))) %>%
  filter(sex=="F", group=="MDD.BPD")

# both/ either L-A or L-R sig
sig.both = left_join(sig.lr, mutate(sig.la[,c("gene_id","gene_name","sex","group","dir")], adj_sig=TRUE), 
                     by=c("gene_id","gene_name","sex","group","dir")) %>%
  mutate(adj_sig=ifelse(is.na(adj_sig), F, T))

sig.both = bind_rows(sig.both, 
          filter(sig.la, adj.P.Val<.01, !gene_id %in% sig.both$gene_id) %>%
            mutate(cluster= "L-A only", seurat_label_f= "L-A only", adj_sig=T)) %>%
  mutate(fill_color= factor(paste(adj_sig, dir), levels=c("TRUE decreased","TRUE increased",
                                                          "FALSE decreased","FALSE increased"),
                            labels=c("Layer-adj. down","Layer-adj. up",
                                     "Layer-restr. down","Layer-restr. up")),
         seurat_label_f=factor(cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb", "L-A only"),
                               labels=c("M.V","Astro","L2.3","L4","L5","L6","Oligo","Inhb", "L-A"))
  )


if (!"MDD-BPD-F_label-transfer-age_LA-LR-sig-DEGs.csv" %in% list.files("processed-data/07_dx_DE/")) {
  write.csv(sig.both, "processed-data/07_dx_DE/MDD-BPD-F_label-transfer-age_LA-LR-sig-DEGs.csv", row.names=F)
}

### volcano plot of L-R and L-A results
plist = list()

plot.df = left_join(filter(restr.results, group=="MDD.BPD", sex=="F"),
                   sig.both[,c("gene_name","gene_id","group","sex","seurat_label_f","dir","adj_sig")]) %>%
  mutate(adj_sig= ifelse(is.na(adj_sig), F, adj_sig))
plist[[1]] = ggplot(plot.df, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3, color="grey")+
  geom_point(data=filter(plot.df, !is.na(dir)), size=.3, color="black")+
  geom_point(data=filter(plot.df, adj_sig==T), size=.3, color="red2")+
  facet_wrap(vars(seurat_label_f), ncol=4)+
  xlim(-ceiling(max(plot.df$logFC)), ceiling(max(plot.df$logFC)))+
  labs(title="Red= L-R and L-A sig.; Black= L-R sig. only")+
  theme_bw()

### bar plot of L-R and L-A results
fill_pal = c("Layer-adj. down"="skyblue","Layer-adj. up"="tomato",
             "Layer-restr. down"="skyblue3", "Layer-restr. up"="red3")

plot.df = group_by(sig.both, seurat_label_f, dir, fill_color, .drop=F) %>% tally() %>%
  mutate(cluster= factor(seurat_label_f, levels=rev(levels(sig.both$seurat_label_f))))

plist[[2]] = ggplot(plot.df, aes(y=cluster, x=n, fill=fill_color))+
  geom_bar(data=filter(plot.df, dir=="decreased", fill_color %in% c("Layer-adj. down","Layer-restr. down")), 
	aes(x=-n), stat="identity", 
           position= position_dodge(preserve="single"), color="black", linewidth=.3)+
  geom_bar(data=filter(plot.df, dir=="increased", fill_color %in% c("Layer-adj. up","Layer-restr. up")), stat="identity", 
           position= position_dodge(preserve="single"), color="black", linewidth=.3)+
  #facet_grid(rows=vars(group), cols=vars(sex))+
  scale_fill_manual(values=fill_pal)+
  xlim(-max(plot.df$n), max(plot.df$n))+
  labs(title="L-A in y-axis= L-A sig. only (no L-R sig.)", 
       fill="",
       x="# sig. genes (sign indicates directionality)", y="")+
  theme_bw()

#either direction
restr.name = unique(filter(sig.both, #dir=="decreased", 
                           cluster!="L-A only")$gene_name)
adj.name = unique(filter(sig.both, #dir=="decreased", 
                         adj_sig==T)$gene_name)
sig.name= union(restr.name, adj.name)

#consistent with NTC results
incr.bpd = unique(filter(sig.both, dir=="increased")$gene_name)
incr.mdd = unique(filter(sig.both, dir=="decreased")$gene_name)

ntc.mdd.dn = bind_rows(filter(restr.results, cluster!="Oligo", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="M", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="F", adj.P.Val<.0001)) %>%
  mutate(dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased","increased"))) %>%
  filter(sex=="F", group=="NTC.MDD", dir=="decreased") %>%
  bind_rows(filter(adj.results, sex=="F", group=="NTC.MDD", adj.P.Val<.01, dir=="decreased"))
decr.mdd.consistent = intersect(ntc.mdd.dn$gene_name, incr.bpd)


ntc.mdd.up = bind_rows(filter(restr.results, cluster!="Oligo", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="M", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="F", adj.P.Val<.0001)) %>%
  mutate(dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased","increased"))) %>%
  filter(sex=="F", group=="NTC.MDD", dir=="increased") %>%
  bind_rows(filter(adj.results, sex=="F", group=="NTC.MDD", adj.P.Val<.01, dir=="increased"))
incr.mdd.consistent = intersect(ntc.mdd.up$gene_name, incr.mdd)


ntc.bpd.up = bind_rows(filter(restr.results, cluster!="Oligo", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="M", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="F", adj.P.Val<.0001)) %>%
  mutate(dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased","increased"))) %>%
  filter(sex=="F", group=="NTC.BPD", dir=="increased") %>%
  bind_rows(filter(adj.results, sex=="F", group=="NTC.BPD", adj.P.Val<.01, dir=="increased"))
incr.bpd.consistent = intersect(ntc.bpd.up$gene_name, incr.bpd)


ntc.bpd.dn = bind_rows(filter(restr.results, cluster!="Oligo", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="M", adj.P.Val<.01),
                       filter(restr.results, cluster=="Oligo", sex=="F", adj.P.Val<.0001)) %>%
  mutate(dir= factor(sign(logFC), levels=c(-1,1), labels=c("decreased","increased"))) %>%
  filter(sex=="F", group=="NTC.BPD", dir=="decreased") %>%
  bind_rows(filter(adj.results, sex=="F", group=="NTC.BPD", adj.P.Val<.01, dir=="decreased"))
decr.bpd.consistent = intersect(ntc.bpd.dn$gene_name, incr.mdd)

### BPD up first because only 5 genes
# BPD up genes summary table 
cat("\n\nMDD.BPD F increased...\n")
restr.name1 = unique(filter(sig.both, dir=="decreased", cluster!="L-A only")$gene_name)
adj.name1 = unique(filter(sig.both, dir=="decreased", adj_sig==T)$gene_name)


restr.name2 = unique(filter(sig.both, dir=="increased", cluster!="L-A only")$gene_name)
adj.name2 = unique(filter(sig.both, dir=="increased", adj_sig==T)$gene_name)


gt = as.data.frame(list(group=rep("MDD.BPD", 8), 
                        sex=rep("F", 8),
                        dir=c("decreased","decreased","decreased","decreased",
                              "increased","increased","increased","increased"),
                        model=c("L-A","L-R","overlap","total",
                                "L-A","L-R","overlap","total"),
                        n_sig=c(length(adj.name1), length(restr.name1), length(intersect(restr.name1, adj.name1)), length(union(restr.name1, adj.name1)),
                                length(adj.name2), length(restr.name2), length(intersect(restr.name2, adj.name2)), length(union(restr.name2, adj.name2))
                                )
                        )
                   )

t1 <- ttheme_default(core=list(
        bg_params = list(fill=as.character(factor(gt$dir, levels=c("decreased", "increased"), labels=c("grey90", "tomato")))
			)
        ))

plist[[3]] = tableGrob(gt, rows = NULL, theme=t1)
gt

### L-R genes summary table
gt = group_by(sig.both, group, sex, seurat_label_f, dir, adj_sig) %>% tally() %>%
  tidyr::pivot_wider(names_from="adj_sig", values_from="n", values_fill=0, names_prefix = "adj_sig_")

t2 <- ttheme_default(core=list(
        bg_params = list(fill=as.character(factor(gt$dir, levels=c("decreased", "increased"), labels=c("grey90", "tomato")))
			)
        ))

plist[[4]] = tableGrob(gt, rows = NULL, theme=t2)


#heatmaps
### expression heatmap
row.annot = cbind.data.frame("MDD.BPD L-A sig."= as.character(sig.name %in% adj.name),
                             "MDD.BPD L-R sig."= as.character(sig.name %in% restr.name))
rownames(row.annot) = sig.name
row.annot[incr.bpd, "MDD.BPD L-A sig."] = ifelse(row.annot[incr.bpd,"MDD.BPD L-A sig."]==T, "BPD up", FALSE)
row.annot[incr.bpd, "MDD.BPD L-R sig."] = ifelse(row.annot[incr.bpd,"MDD.BPD L-R sig."]==T, "BPD up", FALSE)


row.annot[incr.mdd, "MDD.BPD L-A sig."] = ifelse(row.annot[incr.mdd,"MDD.BPD L-A sig."]==T, "MDD up", FALSE)
row.annot[incr.mdd, "MDD.BPD L-R sig."] = ifelse(row.annot[incr.mdd,"MDD.BPD L-R sig."]==T, "MDD up", FALSE)

row.annot$NTC.BPD = "FALSE"
row.annot[incr.bpd.consistent, "NTC.BPD"] = "BPD up"
row.annot[decr.bpd.consistent, "NTC.BPD"] = "BPD dn"


annot_colors = list("MDD.BPD L-A sig."=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "BPD up"=cpList$dx.pal[["BPD"]]),
                    "MDD.BPD L-R sig."=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "BPD up"=cpList$dx.pal[["BPD"]]),
                    #"NTC.MDD"=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "MDD dn"=cpList$dx.pal[["NTC"]]),
                    "NTC.BPD"=c("FALSE"="white","BPD up"=cpList$dx.pal[["BPD"]], "BPD dn"=cpList$dx.pal[["NTC"]]))


spe_pseudo$seurat_label2 = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Oligo","Astro",
                                                                    "L2.3","L4","L5","L6","Inhb"),
                                  labels=c("M.V","Oligo","Astro",
                                           "L2.3","L4","L5","L6","Inhb"))

phm = plotGroupedHeatmap(spe_pseudo, features=incr.bpd, swap_rownames="gene_name",
                         group="seurat_label2", cluster_cols=F, center=T, angle_col=45,
                         fontsize_row=8, 
                         annotation_row=row.annot, annotation_colors=annot_colors,
                         main="MDD.BPD F: Increased (scaled logcounts)")

plist[[5]] = phm[[4]]

### L-R t stat heatmap
t1 = filter(restr.results, group=="MDD.BPD", sex=="F", gene_name %in% sig.name) %>% 
  select(gene_name, seurat_label_f, t) %>%
  tidyr::pivot_wider(names_from="seurat_label_f", values_from="t", values_fill=NA)


m1 = as.matrix(t1[,-1])
rownames(m1) = t1$gene_name
v = ceiling(max(abs(m1)))

phm = pheatmap(m1[incr.bpd, c("M.V","Oligo","Astro","L2.3","L4","L5","L6","Inhb")], 
               cluster_rows=T, cluster_cols=F, 
               annotation_row=row.annot, annotation_colors=annot_colors,
               fontsize_row=8, 
               angle_col=45, 
               colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))((v*2)+1),
               breaks=seq(-v, v, length.out=((v*2)+2)),
               legend_breaks = seq(-v, v, by=2),
               main="MDD.BPD F: Increased (t statistics)")

plist[[6]] = phm[[4]]


top.genes = c("TTYH1","TMEM132A","CRTC1", "KIF5A")

# ORA (MDD > BPD) --------------------------------------------------------
# in order to take advantage of the many gmt objects for enrichR but also the easy plotting with clusterProfiler
# i made a custom workflow that pulls the enrichR gmt and formats it for use with clusterProfiler::enricher
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
#dbs[grep("GO",dbs$libraryName),]



### modify enrichR .read_gmt to return it gmt file of enrichr database as a gson 
#https://github.com/wjawaid/enrichR/blob/master/R/functions.R#L232

#### modify this function so it checks a jhpce dir for the gmt file and if not there finds it online and saves it to jhpce
.read_gmt <- function(db) {
  dbs <- listEnrichrDbs()
  if(!db %in% dbs$libraryName) stop("Requested database not in enrichR. Use enrichR::listEnrichrDbs() to check for available resources.")
  gmtDir = "code/enrichR_gmts"
  if(paste0(db,".rda") %in% list.files(gmtDir)) {
    gmt = readRDS(paste0(gmtDir, "/", db,".rda"))
  } else {
    base.address <- getOption("enrichR.base.address")
    url <- paste0(base.address, "geneSetLibrary?mode=text&libraryName=", db)
    tf <- tempfile(pattern = db, fileext = ".gmt")
    cat("   - Download GMT file...\n")
    tryCatch(download.file(url, tf, mode = "w", quiet = TRUE),
             warning = function(warn) { message(warn); message("") },
             error = function(err) { message(err); message("") })
    gmt = read.gmt(tf)
    attr(gmt, 'snapshot') = c("retrieved"=format(Sys.time()), enrichR_version=packageVersion("enrichR"))
    attr(gmt, 'dbs_details') = as.list(dbs[grep(db, dbs$libraryName),])
    saveRDS(gmt, paste0(gmtDir, "/", db, ".rda"))
    cat("   - Saved to:", paste0(gmtDir, "/", db, ".rda"),"\n")
  }

  return(gmt)
}


### GO BP
cat("\n\n>>> GO (BP)...\n")
go.gmt = .read_gmt("GO_Biological_Process_2025")
#separate term and GOID for TERM2NAME
go.term = do.call(rbind.data.frame, strsplit(as.character(go.gmt$term), split=" \\(GO:"))
colnames(go.term) <- c("Term","ID")
go.term$ID = gsub("\\)", "", paste0("GO:",go.term$ID))
go.id.long = go.term$ID
go.term = distinct(go.term[,c("ID","Term")])

#dframe of goID and gene name for TERM2GENE
go.gene = cbind.data.frame("ID"=go.id.long, "geneID"=go.gmt$gene)

go.results = enricher(incr.bpd, #universe=unique(go.gmt$gene), 
                      TERM2GENE = go.gene,
                      TERM2NAME = go.term)
set.seed(123) #reset seed
nrow(filter(go.results@result, p.adjust<.05)) 
go.results@result <- go.results@result[go.results@result$p.adjust<.05, ]


ora.list = list("GO.BP"=go.results@result)


### Reactome
cat("\n\n>>> Reactome...\n")
react.gmt = .read_gmt("Reactome_2022")
#separate term and GOID for TERM2NAME
react.term = do.call(rbind.data.frame, strsplit(as.character(react.gmt$term), split=" R-HSA"))
colnames(react.term) <- c("Term","ID")
react.term$ID = gsub("\\)", "", paste0("R-HSA",react.term$ID))
react.id.long = react.term$ID
react.term = distinct(react.term[,c("ID","Term")])
#dframe of goID and gene name for TERM2GENE
react.gene = cbind.data.frame("ID"=react.id.long, "geneID"=react.gmt$gene)
 
react.results = enricher(incr.bpd, #universe=unique(react.gmt$gene), 
                         TERM2GENE = react.gene,
                         TERM2NAME = react.term)
set.seed(123) #reset seed
nrow(filter(react.results@result, p.adjust<.05))
react.results@result = react.results@result[react.results@result$p.adjust<.05,]

ora.list[["Reactome"]] = react.results@result


### WikiPathways
cat("\n\n>>> WikiPathways...\n")
wiki.gmt = .read_gmt("WikiPathways_2024_Human")
#separate term and GOID for TERM2NAME
wiki.term = do.call(rbind.data.frame, strsplit(as.character(wiki.gmt$term), 
                                               split=" WP"))
colnames(wiki.term) <- c("Term","ID")
wiki.term$ID = paste0("WP",wiki.term$ID)
wiki.id.long = wiki.term$ID
wiki.term = distinct(wiki.term[,c("ID","Term")])
#dframe of goID and gene name for TERM2GENE
wiki.gene = cbind.data.frame("ID"=wiki.id.long, "geneID"=wiki.gmt$gene)

wiki.results = enricher(incr.bpd, #universe=rowData(spe_pseudo)$gene_name, 
                         TERM2GENE = wiki.gene,
                         TERM2NAME = wiki.term)
set.seed(123) #reset seed
nrow(filter(wiki.results@result, p.adjust<.05))
wiki.results@result = wiki.results@result[wiki.results@result$p.adjust<.05,]

ora.list[["WikiPathways"]] = wiki.results@result


### PPI
cat("\n\n>>> PPI...\n")
ppi.gmt = .read_gmt("PPI_Hub_Proteins")
ppi.results = enricher(incr.bpd, #universe=rowData(spe_pseudo)$gene_name, 
                       TERM2GENE = ppi.gmt)
set.seed(123) #reset seed
nrow(filter(ppi.results@result, p.adjust<.05))
ppi.results@result = ppi.results@result[ppi.results@result$p.adjust<.05,]

ora.list[["PPI"]] = ppi.results@result

saveRDS(ora.list, "processed-data/07_dx_DE/MDD-BPD-F_label-transfer-age_ORA-results_increased.rda")

#plot expression

spe_sub = spe_pseudo[,spe_pseudo$sex=="F"]

#sce for cell type comparison
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata")

sce_pseudo$seurat_low.res2 = factor(sce_pseudo$seurat_low.res, 
                                    levels=c("Micro.Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"),
                                    labels=c("M.V","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"))


p1 <- plotExpression(sce_pseudo, features=intersect(top.genes, rowData(sce_pseudo)$gene_name),
                     swap_rownames="gene_name", x="seurat_low.res2", 
                     colour_by="seurat_low.res")+
  scale_color_manual(values=cpList$low.res.bright, guide="none")+
  facet_wrap(vars(Feature), ncol=2)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  labs(title="SZBDMulti-seq snRNA-seq expression")+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())

p2 <- plotExpression(spe_pseudo, features=top.genes,
                     swap_rownames="gene_name", x="seurat_label2", 
                     colour_by="seurat_label")+
  scale_color_manual(values=cpList$transfer.bright, guide="none")+
  facet_wrap(vars(Feature), ncol=2)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  labs(title="MBv SRT expression")+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())


plot.genes = top.genes
for (j in plot.genes) {
  colData(spe_sub)[[j]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==j,]
}

plot.df = as.data.frame(colData(spe_sub)[,c("condition","seurat_label","seurat_label2", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

df = filter(adj.results, gene_name %in% plot.genes, group=="MDD.BPD", sex=="F") %>%
  select(group, sex, gene_name, adj.P.Val) %>%
  mutate(adj.P.Val= format(adj.P.Val, scientific=T, digits=2),
         gene_label = paste0(gene_name,"\nL-A adj. p=\n", adj.P.Val))

f_labels = c(df$gene_label, levels(plot.df$seurat_label2))
names(f_labels) = c(df$gene_name, levels(plot.df$seurat_label2))



text.df1 = group_by(plot.df, key_genes) %>% slice_min(n=1, logcounts, with_ties = F) #%>%
text.df2 = left_join(filter(sig.lr, gene_name %in% plot.genes) %>%
                       select(gene_name, seurat_label_f, adj.P.Val),
                     text.df1[,c("key_genes","logcounts")],
                     by=c("gene_name"="key_genes")) %>%
  mutate(adj.P.Val= paste0("L-R adj.p=\n",format(adj.P.Val, scientific=T, digits=2)),
         gene_name= factor(gene_name, levels=plot.genes))

colnames(text.df2) = c("key_genes","seurat_label2","adj.P.Val","logcounts")


p3 = ggplot(plot.df, aes(x=condition, y=logcounts, color=condition))+
  ggbeeswarm::geom_quasirandom(width=.4)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  geom_text(data= text.df2, aes(x=1, label=adj.P.Val), color="black", size=3, hjust=0, vjust=0)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+
  facet_grid(cols=vars(seurat_label2), rows=vars(key_genes), scales="free_y",
             labeller = as_labeller(f_labels))+
  labs(x="", title="MDD.BPD F Increased: Select genes", subtitle="Only female samples plotted")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y.right = element_text(angle=0, face="italic"))


plist[[7]] = arrangeGrob(p1, p2, rasterize(p3, dpi=200), layout_matrix=rbind(c(1,2), c(3,3), c(3,3)), top="BPD > MDD DEGs")


### MDD up now
cat("\n\nMDD.BPD F decreased...\n")
# MDD up genes summary table 
restr.name1 = unique(filter(sig.both, dir=="decreased", cluster!="L-A only")$gene_name)
adj.name1 = unique(filter(sig.both, dir=="decreased", adj_sig==T)$gene_name)


restr.name2 = unique(filter(sig.both, dir=="increased", cluster!="L-A only")$gene_name)
adj.name2 = unique(filter(sig.both, dir=="increased", adj_sig==T)$gene_name)


gt = as.data.frame(list(group=rep("MDD.BPD", 8),
                        sex=rep("F", 8),
                        dir=c("decreased","decreased","decreased","decreased",
                              "increased","increased","increased","increased"),
                        model=c("L-A","L-R","overlap","total",
                                "L-A","L-R","overlap","total"),
                        n_sig=c(length(adj.name1), length(restr.name1), length(intersect(restr.name1, adj.name1)), length(union(restr.name1, adj.name1)),
                                length(adj.name2), length(restr.name2), length(intersect(restr.name2, adj.name2)), length(union(restr.name2, adj.name2))
                                )
                        )
                   )

t1 <- ttheme_default(core=list(
        bg_params = list(fill=as.character(factor(gt$dir, levels=c("decreased", "increased"), labels=c("dodgerblue", "grey90")))
                        )
        ))

plist[[8]] = tableGrob(gt, rows = NULL, theme=t1)
gt

### L-R genes summary table
gt = group_by(sig.both, group, sex, seurat_label_f, dir, adj_sig) %>% tally() %>%
  tidyr::pivot_wider(names_from="adj_sig", values_from="n", values_fill=0, names_prefix = "adj_sig_")

t2 <- ttheme_default(core=list(
        bg_params = list(fill=as.character(factor(gt$dir, levels=c("decreased", "increased"), labels=c("dodgerblue", "grey90")))
                        )
        ))

plist[[9]] = tableGrob(gt, rows = NULL, theme=t2)


#heatmaps
### expression heatmap
row.annot = cbind.data.frame("MDD.BPD L-A sig."= as.character(sig.name %in% adj.name),
                             "MDD.BPD L-R sig."= as.character(sig.name %in% restr.name))
rownames(row.annot) = sig.name
row.annot[incr.bpd, "MDD.BPD L-A sig."] = ifelse(row.annot[incr.bpd,"MDD.BPD L-A sig."]==T, "BPD up", FALSE)
row.annot[incr.bpd, "MDD.BPD L-R sig."] = ifelse(row.annot[incr.bpd,"MDD.BPD L-R sig."]==T, "BPD up", FALSE)


row.annot[incr.mdd, "MDD.BPD L-A sig."] = ifelse(row.annot[incr.mdd,"MDD.BPD L-A sig."]==T, "MDD up", FALSE)
row.annot[incr.mdd, "MDD.BPD L-R sig."] = ifelse(row.annot[incr.mdd,"MDD.BPD L-R sig."]==T, "MDD up", FALSE)

row.annot$NTC.BPD <- NULL
row.annot$NTC.MDD = "FALSE"
row.annot[incr.mdd.consistent, "NTC.MDD"] = "MDD up"
row.annot[decr.mdd.consistent, "NTC.MDD"] = "MDD dn"

annot_colors = list("MDD.BPD L-A sig."=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "BPD up"=cpList$dx.pal[["BPD"]]),
                    "MDD.BPD L-R sig."=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "BPD up"=cpList$dx.pal[["BPD"]]),
                    "NTC.MDD"=c("FALSE"="white","MDD up"=cpList$dx.pal[["MDD"]], "MDD dn"=cpList$dx.pal[["NTC"]])
                    #"NTC.BPD"=c("FALSE"="white","BPD up"=cpList$dx.pal[["BPD"]], "BPD dn"=cpList$dx.pal[["NTC"]])
                    )

phm = plotGroupedHeatmap(spe_pseudo, features=incr.mdd, swap_rownames="gene_name",
                         group="seurat_label2", cluster_cols=F, center=T, angle_col=45,
                         fontsize_row=6, 
                         annotation_row=row.annot, annotation_colors=annot_colors,
                         main="MDD.BPD F: Decreased (scaled logcounts)")

plist[[10]] = phm[[4]]

### L-R t stat heatmap
phm = pheatmap(m1[incr.mdd, c("M.V","Oligo","Astro","L2.3","L4","L5","L6","Inhb")], 
               cluster_rows=T, cluster_cols=F, 
               annotation_row=row.annot, annotation_colors=annot_colors,
               fontsize_row=6, 
               angle_col=45, 
               colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))((v*2)+1),
               breaks=seq(-v, v, length.out=((v*2)+2)),
               legend_breaks = seq(-v, v, by=2),
               main="MDD.BPD F: Decreased (t statistics)")


plist[[11]] = phm[[4]]


# ORA (MDD > BPD) --------------------------------------------------------

top.genes = intersect(union(adj.name, restr.name), incr.mdd.consistent)
cat("\nTop genes are genes that are MDD > BPD L-A or L-R and also MDD > NTC L-R:\n")
top.genes

### GO BP
cat("\n\n>>> GO (BP)...\n")
go.gmt = .read_gmt("GO_Biological_Process_2025")
#separate term and GOID for TERM2NAME
go.term = do.call(rbind.data.frame, strsplit(as.character(go.gmt$term), split=" \\(GO:"))
colnames(go.term) <- c("Term","ID")
go.term$ID = gsub("\\)", "", paste0("GO:",go.term$ID))
go.id.long = go.term$ID
go.term = distinct(go.term[,c("ID","Term")])

#dframe of goID and gene name for TERM2GENE
go.gene = cbind.data.frame("ID"=go.id.long, "geneID"=go.gmt$gene)

go.results = enricher(incr.mdd, #universe=unique(go.gmt$gene),
                      TERM2GENE = go.gene,
                      TERM2NAME = go.term)
set.seed(123) #reset seed
nrow(filter(go.results@result, p.adjust<.05))
go.results@result <- go.results@result[go.results@result$p.adjust<.05, ]


ora.list = list("GO.BP"=go.results@result)


cat("\n>>>>> Genes repeated in most GO BP terms:\n")
tidyr::separate_rows(go.results@result, geneID, sep="/") %>% 
  group_by(geneID) %>% tally() %>% filter(n>10) %>% arrange(desc(n))

go.results <- pairwise_termsim(go.results)
set.seed(123) #reset seed
baseplot = emapplot(go.results, showCategory=nrow(go.results@result))+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed (not necessary here)


modify_cnetplot <- function(cnet_plot, term_list, gene_list=NULL) {
  m = ggplot_build(cnet_plot)
  
  m$data[[4]]$colour = ifelse(m$data[[4]]$label %in% top.genes, "red4", "black")
  m$data[[4]]$fontface = ifelse(m$data[[4]]$label %in%  adj.name, "bold.italic", "italic")
  m$data[[4]][1:length(term_list), "fontface"] = "plain"
  m$data[[4]][1:length(term_list), "size"] = 2
  #wrap category text
  m$data[[4]]$label = stringr::str_wrap(m$data[[4]]$label, width = 30)
  
  #highlight new genes for reactome and wiki cnetplots
  if(!is.null(gene_list)) {
    tmp = m$data[[4]]
    tmp$color = ifelse(tmp$label %in% union(term_list, gene_list), "keep", "white")
    tmp2 = merge(m$data[[3]], tmp[,c("x","y","color")], sort=F) #nodes
    m$data[[3]]$colour <- ifelse(tmp2$color=="white", "grey85", tmp2$colour)
    m$data[[3]]$size <- ifelse(tmp2$color=="white", 1, tmp2$size)
    m$data[[4]][!m$data[[4]]$label %in% gene_list, "size"] = 2
  }

  return(ggplot_gtable(m))
}

modify_emmapplot <- function(emmap_plot, term_list) {
  m = ggplot_build(emmap_plot)
  tmp = m$data[[3]]
  tmp$color = ifelse(tmp$label %in% term_list, "keep", "white")
  #modify nodes
  tmp2 = merge(m$data[[2]], tmp[,c("x","y","color")], sort=F) #nodes
  m$data[[2]]$size <- ifelse(tmp2$color=="white", 1, tmp2$size)
  m$data[[2]]$colour <- ifelse(tmp2$color=="white", "grey85", tmp2$colour)
  #modify edges
  tmp3 = merge(m$data[[1]], tmp[,c("x","y","color")], sort=F) #edges
  m$data[[1]]$colour <- ifelse(tmp3$color=="white", "grey85", tmp3$colour)
  m$data[[1]]$linewidth <- ifelse(tmp3$color=="white", .2, tmp3$linewidth)
  #modify text
  m$data[[3]]$size = 2
  m$data[[3]]$label = ifelse(!m$data[[3]]$label %in% term_list, "", m$data[[3]]$label)
  m$data[[3]]$label = stringr::str_wrap(m$data[[3]]$label, width = 30)
  
  return(ggplot_gtable(m))
}

#focus on s100a9 terms
s100a9.terms = tidyr::separate_rows(go.results@result, geneID, sep="/") %>% 
  filter(geneID=="S100A9") %>% pull(Description) %>% unique()
go.s100 = go.results
go.s100@result = go.s100@result[go.s100@result$Description %in% s100a9.terms, ]

p = cnetplot(go.s100, showCategory=nrow(go.s100@result), layout="fr")+
  labs(title="MDD.BPD F: Decreased", subtitle="GO (BP) results with S100A9 (n= 24 sig. terms)\nRed= top, consistent genes; Bold= L-A sig.genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed

plist[[12]] = modify_cnetplot(p, s100a9.terms)


plist[[13]] = modify_emmapplot(baseplot+labs(title="MDD.BPD F: Decreased", subtitle="GO (BP) results with S100A9 (n= 24 sig. terms)"),
                              s100a9.terms)


#focus on CD74
cd74.terms = tidyr::separate_rows(go.results@result, geneID, sep="/") %>% 
  filter(geneID=="CD74") %>% pull(Description) %>% unique()
go.cd = go.results
go.cd@result = go.cd@result[go.cd@result$Description %in% cd74.terms, ]

p = cnetplot(go.cd, showCategory=nrow(go.cd@result), layout="fr")+
  labs(title="MDD.BPD F: Decreased", subtitle="GO (BP) results with CD74 (n= 20 sig. terms)\nRed= top, consistent genes; Bold= L-A sig. genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed

plist[[14]] = modify_cnetplot(p, cd74.terms)

plist[[15]] = modify_emmapplot(baseplot+labs(title="MDD.BPD F: Decreased", 
                                             subtitle="GO (BP) results with CD74 (n= 20 sig. terms)"),
                               cd74.terms)


#other terms 
cat("\n\n>>>>> GP BP terms without S100A9 or CD74:\n")
go.other = go.results
go.other@result = go.other@result[!go.other@result$Description %in% union(s100a9.terms, cd74.terms), ]
nrow(go.other@result)

p = cnetplot(go.other, showCategory=nrow(go.other@result), layout="fr",
             size_category=.5, cex_label_category=.5)+
  labs(title="MDD.BPD F: Decreased", subtitle="GO (BP) results without S100A9 or CD74 (n= 17 sig. terms)\nRed= top, consistent genes; Bold= L-A sig. genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed 

plist[[16]] = modify_cnetplot(p, go.other@result$Description)

plist[[17]] = modify_emmapplot(baseplot+labs(title="MDD.BPD F: Decreased", 
                                             subtitle="GO (BP) results without S100A9 or CD74 (n= 17 sig. terms)"), go.other@result$Description)




#other gmts
cat("\n\n>>> Reactome...\n")
react.gmt = .read_gmt("Reactome_2022")
#separate term and GOID for TERM2NAME
react.term = do.call(rbind.data.frame, strsplit(as.character(react.gmt$term), split=" R-HSA"))
colnames(react.term) <- c("Term","ID")
react.term$ID = gsub("\\)", "", paste0("R-HSA",react.term$ID))
react.id.long = react.term$ID
react.term = distinct(react.term[,c("ID","Term")])
#dframe of goID and gene name for TERM2GENE
react.gene = cbind.data.frame("ID"=react.id.long, "geneID"=react.gmt$gene)


react.results = enricher(incr.mdd,
                         #adj.name, #universe=unique(react.gmt$gene), 
                         TERM2GENE = react.gene,
                         TERM2NAME = react.term)
set.seed(123) #reset seed
nrow(filter(react.results@result, p.adjust<.05)) 
react.results@result <- react.results@result[react.results@result$p.adjust<.05, ]

ora.list[["Reactome"]] = react.results@result

#plot new only
cat("\n>>>>> Reactome terms with genes not present in GO (BP) results:\n")
tmp1 = tidyr::separate_rows(go.results@result, geneID, sep="/")
tmp = tidyr::separate_rows(react.results@result, geneID, sep="/")
react.genes = setdiff(tmp$geneID, tmp1$geneID)
react.terms = unique(filter(tmp, geneID %in% react.genes)$Description)
length(react.terms)

react.results2 = react.results
react.results2@result = react.results2@result[react.results2@result$Description %in% react.terms, ]
p = cnetplot(react.results2, showCategory=nrow(react.results2@result), layout="fr",
             size_category=.5, cex_label_category=.5)+
  labs(title="MDD.BPD F: Decreased", subtitle="Reactome results with genes not present in GO (BP) results (n= 8 sig. terms)\nRed= top, consistent genes; Bold= L-A sig. genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed

plist[[18]] = modify_cnetplot(p, react.terms, gene_list=react.genes)


#wikipathways
cat("\n\n>>> WikiPathways...\n")
wiki.gmt = .read_gmt("WikiPathways_2024_Human")
#separate term and GOID for TERM2NAME
wiki.term = do.call(rbind.data.frame, strsplit(as.character(wiki.gmt$term), 
                                               split=" WP"))
colnames(wiki.term) <- c("Term","ID")
wiki.term$ID = paste0("WP",wiki.term$ID)
wiki.id.long = wiki.term$ID
wiki.term = distinct(wiki.term[,c("ID","Term")])
#dframe of goID and gene name for TERM2GENE
wiki.gene = cbind.data.frame("ID"=wiki.id.long, "geneID"=wiki.gmt$gene)

wiki.results = enricher(incr.mdd,
                        #adj.name, #universe=rowData(spe_pseudo)$gene_name, 
                        TERM2GENE = wiki.gene,
                        TERM2NAME = wiki.term)
set.seed(123) #reset seed
nrow(filter(wiki.results@result, p.adjust<.05)) 
wiki.results@result <- wiki.results@result[wiki.results@result$p.adjust<.05, ]

ora.list[["WikiPathways"]] = wiki.results@result

#plot only new
cat("\n>>>>> WikiPathways terms with genes not present in GO (BP) results:\n")
tmp = tidyr::separate_rows(wiki.results@result, geneID, sep="/")
wiki.genes = setdiff(tmp$geneID, tmp1$geneID)
wiki.terms = unique(filter(tmp, geneID %in% wiki.genes)$Description)
length(wiki.terms) 

wiki.results2 = wiki.results
wiki.results2@result = wiki.results2@result[wiki.results2@result$Description %in% wiki.terms, ]
p = cnetplot(wiki.results2, showCategory=nrow(wiki.results2@result), layout="fr",
             size_category=.5, cex_label_category=.5)+
  labs(title="MDD.BPD F: Decreased", subtitle="WikiPathways results with genes not present in GO (BP) results (n= 8 sig. terms)\nRed= top, consistent genes; Bold= L-A sig. genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed

plist[[19]] = modify_cnetplot(p, wiki.terms, gene_list=wiki.genes)


#ppi
cat("\n\n>>> PPI...\n")
ppi.gmt = .read_gmt("PPI_Hub_Proteins")
head(ppi.gmt)
ppi.results = enricher(sig.name,
                       #adj.name, #universe=rowData(spe_pseudo)$gene_name, 
                       TERM2GENE = ppi.gmt)
set.seed(123) #reset seed
nrow(filter(ppi.results@result, p.adjust<.05))
ppi.results@result <- ppi.results@result[ppi.results@result$p.adjust<.05, ]

ora.list[["PPI"]] =  ppi.results@result

saveRDS(ora.list, "processed-data/07_dx_DE/MDD-BPD-F_label-transfer-age_ORA-results_decreased.rda")

### key genes expression
cat("\n\nPlot expression of key genes grouped by function...\n")

select.genes = list("Microglia"=c("S100A9","S100A8","FCER1G","CD74"),
	"Microglial signaling"=c("SPP1","RGS1","IFITM2"),
	"Complement Component"=c("C1QB","C1QC","C1QA","C3AR1"),
	"Actin and Rac"=c("RHOG","TPM4","DTNA","AIF1")
)

for (i in names(select.genes)) {

p1 <- plotExpression(sce_pseudo, features=intersect(select.genes[[i]], rowData(sce_pseudo)$gene_name),
               swap_rownames="gene_name", x="seurat_low.res2", 
               colour_by="seurat_low.res")+
  scale_color_manual(values=cpList$low.res.bright, guide="none")+
  facet_wrap(vars(Feature), ncol=2)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  labs(title="SZBDMulti-seq snRNA-seq expression")+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())

p2 <- plotExpression(spe_pseudo, features=select.genes[[i]],
               swap_rownames="gene_name", x="seurat_label2", 
               colour_by="seurat_label")+
  scale_color_manual(values=cpList$transfer.bright, guide="none")+
  facet_wrap(vars(Feature), ncol=2)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  labs(title="MBv SRT expression")+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())


plot.genes = select.genes[[i]]
for (j in plot.genes) {
  colData(spe_sub)[[j]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==j,]
}

plot.df = as.data.frame(colData(spe_sub)[,c("condition","seurat_label","seurat_label2", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

df = filter(adj.results, gene_name %in% plot.genes, group=="MDD.BPD", sex=="F") %>%
  select(group, sex, gene_name, adj.P.Val) %>%
  mutate(adj.P.Val= format(adj.P.Val, scientific=T, digits=2),
         gene_label = paste0(gene_name,"\nL-A adj. p=\n", adj.P.Val))

f_labels = c(df$gene_label, levels(plot.df$seurat_label2))
names(f_labels) = c(df$gene_name, levels(plot.df$seurat_label2))


text.df1 = group_by(plot.df, key_genes) %>% slice_max(n=1, logcounts, with_ties = F) #%>%
text.df2 = left_join(filter(sig.lr, gene_name %in% plot.genes) %>%
                       select(gene_name, seurat_label_f, adj.P.Val),
            text.df1[,c("key_genes","logcounts")],
            by=c("gene_name"="key_genes")) %>%
  mutate(adj.P.Val= paste0("L-R adj.p=\n",format(adj.P.Val, scientific=T, digits=2)),
         gene_name= factor(gene_name, levels=plot.genes))

colnames(text.df2) = c("key_genes","seurat_label2","adj.P.Val","logcounts")

p3 = ggplot(plot.df, aes(x=condition, y=logcounts, color=condition))+
  ggbeeswarm::geom_quasirandom(width=.4)+
  geom_violin(draw_quantiles = c(.5), fill="transparent", color="black", scale="width")+
  geom_text(data= text.df2, aes(x=1, label=adj.P.Val), color="black", size=3, hjust=0, vjust=1)+
  scale_color_manual(values=cpList$dx.pal, guide="none")+
  facet_grid(cols=vars(seurat_label2), rows=vars(key_genes), scales="free_y",
             labeller = as_labeller(f_labels))+
  labs(x="", title=paste("MDD.BPD F (MDD > BPD):", i), subtitle="Only female samples plotted")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y.right = element_text(angle=0, face="italic"))

plist[[length(plist)+1]] = arrangeGrob(p1, p2, rasterize(p3, dpi=200), layout_matrix=rbind(c(1,2), c(3,3), c(3,3)), top=paste(i, "DEGs"))

}


cat("\n\nSave MDD.BPD F report: plots/07_dx_DE/prelim-report_MDD-BPD-F_label-transfer-age_both.pdf\n")

pdf(file="plots/07_dx_DE/prelim-report_MDD-BPD-F_label-transfer-age_both.pdf", width=8.5, height=11)
grid.arrange(rasterize(plist[[1]], layer="point", dpi=200), plist[[2]], ncol=1)
grid.arrange(plist[[3]], plist[[4]], ncol=1)
grid.arrange(plist[[5]], plist[[6]], ncol=1)
plot(plist[[7]])
grid.arrange(plist[[8]], plist[[9]], ncol=1)
grid.arrange(plist[[10]], plist[[11]], ncol=1)
plot(plist[[12]])
plot(plist[[13]])
plot(plist[[14]])
plot(plist[[15]])
plot(plist[[16]])
plot(plist[[17]])
plot(plist[[18]])
plot(plist[[19]])
plot(plist[[20]])
plot(plist[[21]])
plot(plist[[22]])
plot(plist[[23]])
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
