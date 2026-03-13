results_set="smoothed-k9-1663"
la.degs_sm = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs_sm = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
sm.degs = union(la.degs_sm$gene_name, lr.degs_sm$gene_name)
length(sm.degs) #649

sm.degs2 = union(la.degs_sm$gene_name[la.degs_sm$n_ttest_sig>0], lr.degs_sm$gene_name[lr.degs_sm$n_ttest_sig>0])
length(sm.degs2) #387

results_set="seurat-pc30"
la.degs_se = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
lr.degs_se = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                             "_dx-sex_degs-F-test-t-test.csv"))
se.degs = union(la.degs_se$gene_name, lr.degs_se$gene_name)
length(se.degs) #627

se.degs2 = union(la.degs_se$gene_name[la.degs_se$n_ttest_sig>0], lr.degs_se$gene_name[lr.degs_se$n_ttest_sig>0])
length(se.degs2) #408

mbv.degs = sort(union(sm.degs, se.degs))
length(mbv.degs)

fileConn<-file("raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt")
writeLines(mbv.degs, fileConn)
close(fileConn)
