setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
library(mclust)
library(ggplot2)
library(dplyr)
library(gridExtra)
})

res_cons = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res_orig = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)

setdiff(rownames(res_orig), rownames(res_cons))
length(setdiff(rownames(res_cons), rownames(res_orig)))

res_cons$original = "QC removed"
res_cons[rownames(res_orig),"original"] = res_orig$predicted.id
res_cons$original = factor(res_cons$original, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo","QC removed"),
                                       labels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo","QC removed"))
table(res_cons$original, useNA="ifany")

res_cons$original_merged = factor(res_cons$original, levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo","QC removed"),
                           labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo","QC removed"))
table(res_cons$original_merged, useNA="ifany")


res_cons$conservative = factor(res_cons$predicted.id, levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"))
table(res_cons$conservative, useNA="ifany")

res_cons$conservative_merged = factor(res_cons$predicted.id, levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                                  labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))
table(res_cons$conservative_merged, useNA="ifany")

################# pairwise jaccard
jcoef <- function(x, y) { 
  # x is a named T/F vector for reference cluster ID, with names being spot codes
  # y is a named factor, numeric, or character vector of all comparison cluster IDs, with names being spotcodes
  if(class(y)!="factor") y=as.factor(y)
  x.ids = names(x)[x]
  y.list = levels(y)
  names(y.list) = y.list
  #for all levels of comparison factor, return JC with reference
  sapply(y.list, function(z) {
    z.ids = names(y)[y==z]
    length(intersect(x.ids, z.ids))/length(union(x.ids, z.ids))
  })
}

pairwise_jc <- function(source_dataframe, reference_type, compare_type) {
  #source_dataframe is colData; the rownames need to be the obs names
  #reference type and compare type are columns in source_dataframe
  x_values = unique(source_dataframe[[reference_type]])
  #for all values of reference cluster compute JC
  output_list <- lapply(x_values, function(X) {
    test_x = source_dataframe[[reference_type]]==X
    names(test_x) = rownames(source_dataframe)
    test_y = source_dataframe[[compare_type]]
    names(test_y) = rownames(source_dataframe)
    #for reference cluster ID X, compute JC for all comparison clusters
    jc_output <- jcoef(test_x, test_y)
    cbind.data.frame(ref_type=rep(reference_type, length(jc_output)),
                     ref_clus=rep(X, length(jc_output)),
                     comp_type=rep(compare_type, length(jc_output)),
                     comp_clus=names(jc_output),
                     j.coef=as.numeric(jc_output))
  })
  do.call(rbind, output_list)
}


#jaccard plots
jc.df = pairwise_jc(res_cons, "original", "conservative")

p1 <- ggplot(jc.df, aes(x=ref_clus, y=j.coef, label=as.character(comp_clus)))+
  geom_text(size=3)+ylim(0,1)+scale_x_discrete(expand=expansion(mult=.1))+
  labs(title=paste0("vs. Seurat PC30 (original) (ARI= ", 
                    round(adjustedRandIndex(res_cons$original, 
                                            res_cons$conservative), 3),")"), 
       x="Seurat PC30 (original)", y="Jaccard coef.")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=8))


jc.df2 = pairwise_jc(res_cons, "original_merged", "conservative_merged")

p1.s <- ggplot(jc.df2, aes(x=ref_clus, y=j.coef, label=as.character(comp_clus)))+
  geom_text(size=3)+ylim(0,1)+scale_x_discrete(expand=expansion(mult=.1))+
  labs(title=paste0("vs. Seurat PC30 (original, merged) (ARI= ", 
                    round(adjustedRandIndex(res_cons$original_merged, 
                                            res_cons$conservative_merged), 3),")"), 
       x="Seurat PC30 (original, merged)", y="Jaccard coef.")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=8))


#heatmap plots
tmp.df = as.data.frame(table(res_cons[,c("original", "conservative")]))
colnames(tmp.df) = c("reference","query","Freq")
tmp.df = group_by(tmp.df, reference) %>% mutate(Total=sum(Freq)) %>% ungroup() %>%
  mutate(Prop=Freq/Total, query = factor(query, levels=rev(levels(tmp.df$query))))

p2 <- ggplot(tmp.df, aes(y=query, x=reference, fill=Prop))+
  geom_tile(color="grey50", linewidth=.3)+scale_fill_gradient(low="white",high="black", limits=c(0,1))+
  geom_text(data=union(group_by(tmp.df, query) %>% slice_max(n=2, Freq), 
                       group_by(tmp.df, reference) %>% slice_max(n=2, Freq)) %>% 
              filter(Prop>.1) %>%
              mutate(text_value= paste0(round(Freq/1000, 1), "k")), 
            aes(label=text_value), color="red", size=2, fontface="bold")+
  scale_x_discrete(labels=c("Micro/\nVasc", levels(tmp.df$reference)[2:9],"QC\nremoved"))+
  labs(x="Seurat PC30 (original)", y="Seurat PC30 (conservative)", fill="Prop. of\noriginal\ncluster")+
  theme_minimal()+theme(text=element_text(size=8), legend.key.size=unit(6,"pt"),
                        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
                        legend.box.spacing = unit(2,"pt"),
                        legend.margin=margin(0,0,0,0,"pt"),
                        legend.box.margin = margin(0,2,0,2,"pt"))


tmp.df2 = as.data.frame(table(res_cons[,c("original_merged", "conservative_merged")]))
colnames(tmp.df2) = c("reference","query","Freq")
tmp.df2 = group_by(tmp.df2, reference) %>% mutate(Total=sum(Freq)) %>% ungroup() %>%
  mutate(Prop=Freq/Total, query = factor(query, levels=rev(levels(tmp.df2$query))))

p2.s <- ggplot(tmp.df2, aes(y=query, x=reference, fill=Prop))+
  geom_tile(color="grey50", linewidth=.3)+scale_fill_gradient(low="white",high="black", limits=c(0,1))+
  geom_text(data=union(group_by(tmp.df2, query) %>% slice_max(n=2, Freq), 
                       group_by(tmp.df2, reference) %>% slice_max(n=2, Freq)) %>% 
              filter(Prop>.1) %>%
              mutate(text_value= paste0(round(Freq/1000, 1), "k")), 
            aes(label=text_value), color="red", size=2, fontface="bold")+
  scale_x_discrete(labels=c("Micro/\nVasc",levels(tmp.df2$reference)[2:8],"QC\nremoved"))+
  labs(x="Seurat PC30 (original, merged)", y="Seurat PC30 (conservative, merged)", fill="Prop. of\noriginal\ncluster")+
  theme_minimal()+theme(text=element_text(size=8), legend.key.size=unit(6,"pt"),
                        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
                        legend.box.spacing = unit(2,"pt"),
                        legend.margin=margin(0,0,0,0,"pt"),
                        legend.box.margin = margin(0,2,0,2,"pt"))


ggsave("plots/05_clustering/Seurat/compare-clusters_pc30-original-vs-pc30-conservative.pdf", 
       marrangeGrob(list(p1, p2, p1.s, p2.s), nrow=2, ncol=1, top=NULL), height=6, width=4)

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
