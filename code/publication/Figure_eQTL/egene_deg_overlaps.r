library(dplyr)

source("code/09_DEG_GRN/load_DEGs.r")
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

modules = unique(refined.modules$TF)
names(modules) <- modules
modList <- lapply(modules, function(x) c(x, filter(refined.modules, TF==x)$target))


egene.df = read.csv("processed-data/11_eQTL_coloc/seurat/tables/map_significant_pairs.csv.gz")
seurat_labels = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
names(seurat_labels) = seurat_labels


la.df =   mutate(egene.df, cluster= factor(context, levels=rev(seurat_labels))) %>%
  filter(pair_provenance=="cis_and_independent") %>% 
  distinct(cluster, gene_name) %>% 
  left_join(filter(sig.df, source=="se", cluster=="L-A")[,c("gene_name","gene_id","adj.P.Val","sex.group","dir")],
            by=c("gene_name")) %>%
  filter(!is.na(adj.P.Val))

lr.df = mutate(egene.df, cluster= factor(context, levels=rev(seurat_labels))) %>%
  filter(pair_provenance=="cis_and_independent") %>% 
  distinct(cluster, gene_name) %>% 
  left_join(filter(sig.df, source=="se", cluster!="L-A")[,c("gene_name","gene_id","cluster","adj.P.Val","sex.group","dir")],
            by=c("gene_name","cluster")) %>%
  filter(!is.na(adj.P.Val))

strict.de.df = bind_rows(mutate(la.df[,c("cluster","gene_name","gene_id","sex.group")], model="whole-tissue"),
            lr.df[,c("cluster","gene_name","gene_id","sex.group","model")]) %>%
  group_by(cluster, gene_name, gene_id, sex.group) %>% summarise(model=paste(model, collapse=", ")) %>%
  tidyr::pivot_wider(names_from="sex.group", values_from="model", values_fill="ns") %>%
  mutate(cluster=factor(cluster, levels=rev(seurat_labels)))

egene.deg.df = right_join(egene.df[,c("source","pair_provenance","cis_supported","indep_supported","context","gene_id","gene_name",
                       "variant_id","start_distance","ma_samples","ma_count","pval_nominal",
                       "MDD_gwasVar_strict","MDD_gwasVar_exp","MDD_gwasGene",
                       "BD_gwasVar_strict","BD_gwasVar_exp","BD_gwasGene",
                       "SCZD_gwasVar_strict","SCZD_gwasGene")], 
           strict.de.df, by=c("gene_name","gene_id","context"="cluster"))

mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
               "PLP1",
               "GFAP","GLUL","HSPA1A","MT1M",
               "IFITM3","CD74","A2M",
               "GRIN1","CAMK2N1","GAD1")
for(i in mod_subset) {
  egene.deg.df[[i]] = as.numeric(egene.deg.df$gene_name %in% modList[[i]])
}

egene.deg.df


write.csv(egene.deg.df, "processed-data/publication/supp_tables/egene_degs-only.csv", row.names=F)

