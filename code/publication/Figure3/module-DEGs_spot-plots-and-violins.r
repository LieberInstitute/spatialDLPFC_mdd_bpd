suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(dplyr)
  library(ggplot2)
  library(escheR)
  library(ggrastr)
  library(gridExtra)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")

plot.genes = c("CLDN11","SGK1","NEAT1","PGAM2",
               "HSPA1B","GADD45B","MT1X","CEBPD",
               "IFITM3","APOLD1","A2M","ABCG2")

#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
spotdata = spotdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))
spe$smoothed_k9_1663 = factor(spotdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                              labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))

res.pc30 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res.pc30 = res.pc30[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(res.pc30)))
spe$seurat_pc30 = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","Inhb","L4","L5","L6","Oligo"),
                         labels=c("Micro.Vasc","Astro","L2.3","L2.3","Inhb","L4","L5","L6","Oligo"))


spe_sub = spe[,spe$sample_id=="V13B23-308_D1"]

for(i in plot.genes) {
  spe_sub[[i]] = logcounts(spe_sub)[rowData(spe_sub)$gene_name==i,]
}

plist <- lapply(plot.genes, function(x) {
  p = make_escheR(spe_sub) %>%
    add_ground(var="smoothed_k9_1663", point_size=.5, stroke=.3) %>%
    add_fill(var=x, point_size=.5)
  p <- p+scale_color_manual(values=c(cpList$smoothed.light), guide="none")+
    scale_fill_gradient(low="white",high="black")+
    #labs(title=x)+
    theme(#text=element_text(size=10), plot.title=element_text(face="italic"),
      plot.title=element_blank(), legend.position="bottom",
          legend.key.width= unit(14, "pt"), legend.key.height= unit(8, "pt"), 
      legend.text=element_text(size=7),#, margin=margin(0,0,0,0,"pt")),
      legend.title=element_blank(), legend.box.spacing= unit(0,"pt"))
  return(rasterize(p, layers="Point", dpi=300))
})

# violin
results_set="seurat-pc30"
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))

all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright)[c(1:4,8,5:7)]),
                                                                     labels=c("LA","M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")))

#names(plot.genes) <- plot.genes

all.df_filt = filter(all.df, key_genes %in% plot.genes) 
#calc sd/se
all.df_filt2 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% summarise(ypos=mean(logcounts), ysd=sd(logcounts), n=n(), yse=ysd/sqrt(n)) #%>%

#remove outliers
all.df_filt3 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% 
  mutate(iqr=IQR(logcounts), q1=quantile(logcounts, probs=c(.25)), q3= quantile(logcounts, probs=c(.75)),
         out_low=logcounts< (q1-1.5*iqr), out_high= logcounts> (q3+1.5*iqr)) %>%
  filter(out_low==F, out_high==F)


group_by(all.df_filt3, key_genes) %>% summarise(max_log=round(max(logcounts),2))
#ymax1 = ceiling(max(all.df_filt3$logcounts))
ymax1 = 12
vlist = lapply(plot.genes, function(x) {
  tmp = filter(all.df_filt3, key_genes==x)
  tmp2 = filter(all.df_filt2, key_genes==x)
  
  p2 <- ggplot(tmp, aes(x=cluster, y=logcounts))+
    geom_violin(aes(fill=condition), scale="width", color="transparent", position = position_dodge(width=.8), trim=F, bounds=c(0,ymax1))+
    scale_fill_manual(values=cpList$dx.pal, guide="none")+
    facet_grid(cols=vars(sex), labeller= as_labeller(c("F"="Female","M"="Male")))+
    geom_crossbar(data=tmp2, aes(color=condition, y=ypos, ymax=ypos+yse, ymin=ypos-yse), 
                  position = position_dodge(width=.8), linewidth=.3)+
    scale_color_manual(values=c("black","black","black"), guide="none")+
    scale_y_continuous(breaks=c(0,2,4,6,8,10,12))+
    coord_cartesian(ylim=c(0,ymax1))+
    labs(title=x)+
    theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                     text=element_text(size=8), axis.text=element_text(size=6), plot.title=element_text(face="italic"),
                     axis.title.x=element_blank(), axis.title.y=element_blank(),
                     axis.ticks = element_line(linewidth=.2), strip.text.y.left = element_blank(),
                     panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2))
  
  
  return(p2)
})


outlist = c(plist[1:4], vlist[1:4],
            plist[5:8], vlist[5:8],
            plist[9:12], vlist[9:12])

#lay_mat = cbind(1:4,5:8,5:8)
lay_mat = matrix(c(1,5,5,2,6,6,
                   3,7,7,4,8,8), ncol=6, nrow=2, byrow = T)

ggsave(file="plots/publication/Figure3/module-DEGs_spot-plot-and-violin.pdf",
       marrangeGrob(grobs=outlist, layout_matrix=lay_mat, top=NULL),
#       width=5.2, height=8)
	width=10, height=4)
