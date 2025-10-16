setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
library(mclust)
library(ggplot2)
library(dplyr)
library(gridExtra)
})

cdata_cons = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)
cdata_orig = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)

setdiff(rownames(cdata_orig), rownames(cdata_cons))
length(setdiff(rownames(cdata_cons), rownames(cdata_orig)))

cdata_cons$smoothed_k9_1663_f = "QC removed"
cdata_cons[rownames(cdata_orig),"smoothed_k9_1663_f"] = cdata_orig$smoothed_k9_1663_f
cdata_cons$smoothed_k9_1663_f = factor(cdata_cons$smoothed_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI","QC removed"),
                                       labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI","QC removed"))
table(cdata_cons$smoothed_k9_1663_f, useNA="ifany")

cdata_cons$precast_k9_1663_f = "QC removed"
cdata_cons[rownames(cdata_orig),"precast_k9_1663_f"] = cdata_orig$precast_k9_1663_f
cdata_cons$precast_k9_1663_f = factor(cdata_cons$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI","QC removed"),
                                       labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI","QC removed"))
table(cdata_cons$precast_k9_1663_f, useNA="ifany")

cdata_cons = filter(cdata_cons, precast_k7_1626_f!="missing")
cdata_cons$smoothed_k7_1626_f = factor(cdata_cons$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
                         labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))
table(cdata_cons$smoothed_k7_1626_f, useNA="ifany")

cdata_cons$precast_k7_1626_f = factor(cdata_cons$precast_k7_1626_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
                                       labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))
table(cdata_cons$precast_k7_1626_f, useNA="ifany")


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
jc.df = pairwise_jc(cdata_cons, "precast_k9_1663_f", "precast_k7_1626_f")

p1 <- ggplot(jc.df, aes(x=ref_clus, y=j.coef, label=as.character(comp_clus)))+
  geom_text(size=3)+ylim(0,1)+scale_x_discrete(expand=expansion(mult=.1))+
  labs(title=paste0("vs. PRECAST n=1626 k=7 (ARI= ", round(adjustedRandIndex(cdata_cons$precast_k9_1663_f, cdata_cons$precast_k7_1626_f), 3),")"), 
       x="PRECAST n=1663 k=9", y="Jaccard coef.")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=8))


jc.df2 = pairwise_jc(cdata_cons, "smoothed_k9_1663_f", "smoothed_k7_1626_f")

p1.s <- ggplot(jc.df2, aes(x=ref_clus, y=j.coef, label=as.character(comp_clus)))+
  geom_text(size=3)+ylim(0,1)+scale_x_discrete(expand=expansion(mult=.1))+
  labs(title=paste0("vs. PRECAST (smoothed) n=1626 k=7 (ARI= ", 
                    round(adjustedRandIndex(cdata_cons$smoothed_k9_1663_f, 
                                            cdata_cons$smoothed_k7_1626_f), 3),")"), 
       x="PRECAST (smoothed) n=1663 k=9", y="Jaccard coef.")+
  theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1), text=element_text(size=8))


#heatmap plots
tmp.df = as.data.frame(table(cdata_cons[,c("precast_k9_1663_f", "precast_k7_1626_f")]))
colnames(tmp.df) = c("reference","query","Freq")
tmp.df = group_by(tmp.df, reference) %>% mutate(Total=sum(Freq)) %>% ungroup() %>%
  mutate(Prop=Freq/Total, query = factor(query, levels=rev(levels(tmp.df$query))))

p2 <- ggplot(tmp.df, aes(y=query, x=reference, fill=Prop))+
  geom_tile(color="grey50", linewidth=.3)+scale_fill_gradient(low="white",high="black", limits=c(0,1))+
  geom_text(data=#union(group_by(tmp.df, query) %>% slice_max(n=2, Freq), 
                 #      group_by(tmp.df, reference) %>% slice_max(n=2, Freq)) %>% 
              filter(tmp.df, Prop>.1) %>%
              mutate(text_value= paste0(round(Freq/1000, 1), "k")), 
            aes(label=text_value), color="red", size=2, fontface="bold")+
  scale_x_discrete(labels=c(levels(tmp.df$reference)[1:9],"QC\nremoved"))+
  labs(x="PRECAST n=1663 k=9", y="PRECAST n=1626 k=7", fill="Prop. of\nn1663\ncluster")+
  theme_minimal()+theme(text=element_text(size=8), legend.key.size=unit(6,"pt"),
                        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
                        legend.box.spacing = unit(2,"pt"),
                        legend.margin=margin(0,0,0,0,"pt"),
                        legend.box.margin = margin(0,2,0,2,"pt"))


tmp.df2 = as.data.frame(table(cdata_cons[,c("smoothed_k9_1663_f", "smoothed_k7_1626_f")]))
colnames(tmp.df2) = c("reference","query","Freq")
tmp.df2 = group_by(tmp.df2, reference) %>% mutate(Total=sum(Freq)) %>% ungroup() %>%
  mutate(Prop=Freq/Total, query = factor(query, levels=rev(levels(tmp.df2$query))))

p2.s <- ggplot(tmp.df2, aes(y=query, x=reference, fill=Prop))+
  geom_tile(color="grey50", linewidth=.3)+scale_fill_gradient(low="white",high="black", limits=c(0,1))+
  geom_text(data=#union(group_by(tmp.df2, query) %>% slice_max(n=2, Freq), 
                 #      group_by(tmp.df2, reference) %>% slice_max(n=2, Freq)) %>% 
              filter(tmp.df2, Prop>.1) %>%
              mutate(text_value= paste0(round(Freq/1000, 1), "k")), 
            aes(label=text_value), color="red", size=2, fontface="bold")+
  scale_x_discrete(labels=c(levels(tmp.df2$reference)[1:9],"QC\nremoved"))+
  labs(x="PRECAST (smoothed) n=1663 k=9", y="PRECAST n=1626 (smoothed) k=7", fill="Prop. of\nn1663\ncluster")+
  theme_minimal()+theme(text=element_text(size=8), legend.key.size=unit(6,"pt"),
                        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
                        legend.box.spacing = unit(2,"pt"),
                        legend.margin=margin(0,0,0,0,"pt"),
                        legend.box.margin = margin(0,2,0,2,"pt"))


ggsave("plots/05_clustering/PRECAST/compare-clusters_original-vs-conservative.pdf", 
       marrangeGrob(list(p1, p2, p1.s, p2.s), nrow=2, ncol=1, top=NULL), height=6, width=4)

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
