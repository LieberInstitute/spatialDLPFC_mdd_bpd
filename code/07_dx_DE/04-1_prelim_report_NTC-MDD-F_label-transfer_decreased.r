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

sig.la = filter(adj.results, sex=="F", group=="NTC.MDD", adj.P.Val<.01)


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
  filter(sex=="F", group=="NTC.MDD")

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

if (!"NTC-MDD-F_label-transfer-age_LA-LR-sig-DEGs.csv" %in% list.files("processed-data/07_dx_DE/")) {
  write.csv(sig.both, "processed-data/07_dx_DE/NTC-MDD-F_label-transfer-age_LA-LR-sig-DEGs.csv", row.names=F)
}

### genes summary table 
restr.name1 = unique(filter(sig.both, dir=="decreased", cluster!="L-A only")$gene_name)
adj.name1 = unique(filter(sig.both, dir=="decreased", adj_sig==T)$gene_name)


restr.name2 = unique(filter(sig.both, dir=="increased", cluster!="L-A only")$gene_name)
adj.name2 = unique(filter(sig.both, dir=="increased", adj_sig==T)$gene_name)


gt = as.data.frame(list(group=rep("NTC.MDD", 8), 
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

plist = list(tableGrob(gt, rows = NULL, theme=t1))
gt

### L-R genes summary table
gt = group_by(sig.both, group, sex, seurat_label_f, dir, adj_sig) %>% tally() %>%
  tidyr::pivot_wider(names_from="adj_sig", values_from="n", values_fill=0, names_prefix = "adj_sig_")

t2 <- ttheme_default(core=list(
        bg_params = list(fill=as.character(factor(gt$dir, levels=c("decreased", "increased"), labels=c("dodgerblue", "grey90")))
			)
        ))

plist[[2]] = tableGrob(gt, rows = NULL, theme=t2)

### volcano plot of L-R and L-A results
plot.df = left_join(filter(restr.results, group=="NTC.MDD", sex=="F"),
                   sig.both[,c("gene_name","gene_id","group","sex","seurat_label_f","dir","adj_sig")]) %>%
  mutate(adj_sig= ifelse(is.na(adj_sig), F, adj_sig))
plist[[3]] = ggplot(plot.df, aes(x=logFC, y=-log10(adj.P.Val)))+
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

plist[[4]] = ggplot(plot.df, aes(y=cluster, x=n, fill=fill_color))+
  geom_bar(data=filter(plot.df, dir=="decreased", fill_color %in% c("Layer-adj. down","Layer-restr. down")), 
	aes(x=-n), stat="identity", 
           position= position_dodge(preserve="single"), color="black", linewidth=.3)+
  geom_bar(data=filter(plot.df, dir=="increased", fill_color %in% c("Layer-adj. up","Layer-restr. up")), 
	stat="identity", 
           position= position_dodge(preserve="single"), color="black", linewidth=.3)+
  #facet_grid(rows=vars(group), cols=vars(sex))+
  scale_fill_manual(values=fill_pal)+
  xlim(-max(plot.df$n), max(plot.df$n))+
  labs(title="L-A in y-axis= L-A sig. only (no L-R sig.)", 
       fill="",
       x="# sig. genes (sign indicates directionality)", y="")+
  theme_bw()


# NTC.MDD F decreased
cat("\n\nNTC.MDD F decreased...\n")
restr.name = unique(filter(sig.both, dir=="decreased", cluster!="L-A only")$gene_name)
adj.name = unique(filter(sig.both, dir=="decreased", adj_sig==T)$gene_name)
sig.name= union(restr.name, adj.name)

### expression heatmap
row.annot = cbind.data.frame("L-A sig."= as.character(sig.name %in% adj.name),
                             "L-R sig."= as.character(sig.name %in% restr.name))
rownames(row.annot) = sig.name
annot_colors = list("L-A sig."=c("FALSE"="white","TRUE"="black"),
                    "L-R sig."=c("FALSE"="white","TRUE"="black"))

spe_pseudo$seurat_label2 = factor(spe_pseudo$seurat_label, levels=c("Micro.Vasc","Oligo","Astro",
                                                                "L2.3","L4","L5","L6","Inhb"),
                                labels=c("M.V","Oligo","Astro",
                                         "L2.3","L4","L5","L6","Inhb"))
phm = plotGroupedHeatmap(spe_pseudo, features=sig.name, swap_rownames="gene_name",
                                 group="seurat_label2", cluster_cols=F, center=T, angle_col=45,
                                 show_rownames=F, annotation_row=row.annot, annotation_colors=annot_colors,
                         main="NTC.MDD F: Decreased (scaled logcounts)")

plist[[5]] = phm[[4]]


### L-R t stat heatmap
t1 = filter(restr.results, group=="NTC.MDD", sex=="F", gene_name %in% sig.name) %>% 
  select(gene_name, seurat_label_f, t) %>%
  tidyr::pivot_wider(names_from="seurat_label_f", values_from="t", values_fill=NA)


m1 = as.matrix(t1[,-1])
rownames(m1) = t1$gene_name
v = ceiling(max(abs(m1)))
phm = pheatmap(m1[,c("M.V","Oligo","Astro","L2.3","L4","L5","L6","Inhb")], 
         cluster_rows=T, cluster_cols=F, 
         annotation_row=row.annot, annotation_colors=annot_colors,
         show_rownames=F, angle_col=45, 
         colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))((v*2)+1),
         breaks=seq(-v, v, length.out=((v*2)+2)),
         legend_breaks = seq(-v, v, by=2),
         main="NTC.MDD F: Decreased (t statistics)")

plist[[6]] = phm[[4]]


top.genes = phm$tree_row$label[phm$tree_row$order][25:52]



# plot top most consistent genes
phm = pheatmap(m1[top.genes, c("M.V","Oligo","Astro","L2.3","L4","L5","L6","Inhb")], 
               cluster_rows=T, cluster_cols=F, 
               annotation_row=row.annot, annotation_colors=annot_colors,
               #fontsize_row=6, 
               angle_col=45, 
               colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))((v*2)+1),
               breaks=seq(-v, v, length.out=((v*2)+2)),
               legend_breaks = seq(-v, v, by=2), silent=T, 
               main="NTC.MDD F: Top increased genes (t statistics)")

plist[[7]] = phm[[4]]

phm = plotGroupedHeatmap(spe_pseudo, features=top.genes, swap_rownames="gene_name",
                         group="seurat_label2", 
                         cluster_cols=F, cluster_rows=T, center=T, angle_col=45,
                         annotation_row=row.annot, annotation_colors=annot_colors,
                         silent=T, main="NTC.MDD F: Top increased genes (scaled logcounts)")


plist[[8]] = phm[[4]]


#ORA
cat("\n>> ORA...\n")
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()

### modify enrichR .read_gmt to return it gmt file of enrichr database as a gson 
#https://github.com/wjawaid/enrichR/blob/master/R/functions.R#L232
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

go.results = enricher(sig.name, #universe=unique(go.gmt$gene), 
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
 
react.results = enricher(sig.name, #universe=unique(react.gmt$gene), 
                         TERM2GENE = react.gene,
                         TERM2NAME = react.term)
set.seed(123) #reset seed
nrow(filter(react.results@result, p.adjust<.05))
react.results@result <- react.results@result[react.results@result$p.adjust<.05, ]

ora.list[["Reactome"]] = react.results@result


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

wiki.results = enricher(sig.name,
                        #adj.name, #universe=rowData(spe_pseudo)$gene_name, 
                        TERM2GENE = wiki.gene,
                        TERM2NAME = wiki.term)
set.seed(123) #reset seed
nrow(filter(wiki.results@result, p.adjust<.05))
wiki.results@result <- wiki.results@result[wiki.results@result$p.adjust<.05, ]

ora.list[["WikiPathways"]] = wiki.results@result


### merge results and plot
merge_results <- react.results
merge_results@result = rbind(react.results@result, go.results@result, wiki.results@result)
p = cnetplot(merge_results, showCategory=nrow(merge_results@result), layout="fr", 
             size_category=.5, cex_label_category=.5)+
  labs(title="NTC.MDD F: Decreased", subtitle="Green= Reactome; Blue= GO (BP); Pink= WikiPathways\nRed= top sig. genes; Bold= L-A sig. genes")+
  guides("size"=guide_legend("# genes", override.aes = list(color="#B3B3B3")))+
  theme(plot.margin = margin(.5,.5,.5,.5, "cm"))
set.seed(123) #reset seed

m = ggplot_build(p)

#color category nodes by source
tmp = m$data[[4]]
### reactome = palegreen, GO BP = pink
tmp$color = factor(tmp$label, levels=c(react.results@result$Description, go.results@result$Description, wiki.results@result$Description,
                                       setdiff(tmp$label, merge_results@result$Description)),
                   labels=c(rep("palegreen",nrow(react.results@result)), 
                            rep("skyblue",nrow(go.results@result)),
			    rep("pink", nrow(wiki.results@result)), 
                            rep("grey50",nrow(tmp)-nrow(merge_results@result))))
tmp2 = merge(m$data[[2]], tmp[,c("x","y","color")], sort=F)
tmp3 = merge(m$data[[1]], tmp[,c("x","y","color")], sort=F) #edges
m$data[[2]]$colour <- tmp2$color
m$data[[1]]$colour <- tmp3$color #edges

#make gene nodes that are sig. in adjusted a different color
m$data[[4]]$colour = ifelse(m$data[[4]]$label %in% top.genes, "red4", "black")
#make gene nodes italic and key gene nodes bold italic
m$data[[4]]$fontface = ifelse(m$data[[4]]$label %in%  adj.name,
                              "bold.italic", "italic")
m$data[[4]]$fontface = ifelse(m$data[[4]]$label %in% merge_results@result$Description, "plain", m$data[[4]]$fontface)
#make category text smaller
m$data[[4]]$size[1:nrow(merge_results@result)] = 2
#wrap category text
m$data[[4]]$label = stringr::str_wrap(m$data[[4]]$label, width = 30)

plist[[9]] = ggplot_gtable(m)


#ppi
cat("\n\n>>> PPI...\n")
ppi.gmt = .read_gmt("PPI_Hub_Proteins")
ppi.results = enricher(sig.name,
                       #adj.name, #universe=rowData(spe_pseudo)$gene_name, 
                       TERM2GENE = ppi.gmt)
set.seed(123) #reset seed
nrow(filter(ppi.results@result, p.adjust<.05)) #41
ppi.results@result <- ppi.results@result[ppi.results@result$p.adjust<.05, ]

ora.list[["PPI"]] = ppi.results@result


gt = ppi.results@result[, c("Description","GeneRatio","p.adjust","Count","geneID")]
nrow(gt)
gt$geneID = stringr::str_wrap(gsub("/"," ", gt$geneID), width=20)
plist[[10]] = tableGrob(gt, rows = NULL)


saveRDS(ora.list, "processed-data/07_dx_DE/NTC-MDD-F_label-transfer-age_ORA-results_decreased.rda")


### key genes expression
spe_sub = spe_pseudo[,spe_pseudo$sex=="F"]

#sce for cell type comparison
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata")

sce_pseudo$seurat_low.res2 = factor(sce_pseudo$seurat_low.res, 
                                    levels=c("Micro.Vasc","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"),
                                    labels=c("M.V","Astro","Oligo","L2","L3","L4","L5","L6","Inhb"))

select.genes = list("Neuron-neuron signaling"= c("CORT","SST","CRH","VGF"),
	"GABA processing"= c("SLC6A1","GABBR1","GAD2","GABRD"),
	"Postsynaptic signaling"= c("SLC38A5","NCALD","GRM2","NRGN"),
	"Synaptic release"= c("CACNA1B","CPLX1","CPLX2","STX1A","SYNGR1")
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

plot.df = as.data.frame(colData(spe_sub)[,c("condition","seurat_label","seurat_label2",plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

df = filter(adj.results, gene_name %in% plot.genes, group=="NTC.MDD", sex=="F") %>%
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
  labs(x="", title=paste("NTC.MDD F Decreased:", i), subtitle="Only female samples plotted")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y.right = element_text(angle=0, face="italic"))

plist[[length(plist)+1]] = arrangeGrob(p1, p2, rasterize(p3, dpi=200), layout_matrix=rbind(c(1,2), c(3,3), c(3,3)), top=paste(i, "DEGs"))

}


#save ntc.mdd f decreased
cat("\n\nSave NTC.MDD F decreased report: plots/07_dx_DE/prelim-report_NTC-MDD-F_label-transfer-age_decreased.pdf\n")

pdf(file="plots/07_dx_DE/prelim-report_NTC-MDD-F_label-transfer-age_decreased.pdf", width=8.5, height=11)
grid.arrange(plist[[1]], plist[[2]])
grid.arrange(rasterize(plist[[3]], layer="point", dpi=200), plist[[4]], ncol=1)
grid.arrange(plist[[5]], plist[[6]], ncol=2)
grid.arrange(plist[[7]], plist[[8]], ncol=1)
plot(plist[[9]])
grid.arrange(plist[[10]], ncol=1, top="PPI results", bottom="No enriched genes sig. decreased in MDD vs NTC F")
plot(plist[[11]])
plot(plist[[12]])
plot(plist[[13]])
plot(plist[[14]])
dev.off()


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

