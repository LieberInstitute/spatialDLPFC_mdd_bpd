library(dplyr)
library(ggplot2)

### load in F-test sig genes to test subset
ftest = read.table("raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt")

# two mdd options: Jaffe 2022 and Girgenti 2021
gg21 = read.csv("raw-data/bulk_DE_results/Girgenti-2021_ST9.csv")
colnames(gg21)[1:2] = c("gene_id","gene_name")
colnames(gg21)[9:10] = c("logFC","adj.P.Val")
nrow(gg21) # 14996

jf22 = read.csv("raw-data/bulk_DE_results/Jaffe-2022_dlPFC-MDD-only.csv")
colnames(jf22)[c(1:2)] = c("gene_name","gene_id")
colnames(jf22)[6:9] = c("logFC", "t", "P.Value", "adj.P.Val")
nrow(jf22) # 26020

# compare to each other
blk.mdd = inner_join(gg21[,c("gene_id","gene_name","logFC","adj.P.Val")],
	jf22[,c("gene_id","gene_name","logFC","adj.P.Val")],
	by="gene_id", suffix=c("_gg21", "_jf22"))
dim(blk.mdd) # 14414 7
#colSums(is.na(blk.mdd)) # check that no NAs

cor(blk.mdd$logFC_gg21, blk.mdd$logFC_jf22)
# 0.2371127
plot(blk.mdd$logFC_gg21, blk.mdd$logFC_jf22)
# not good looking

plot(-log10(blk.mdd$adj.P.Val_gg21), -log10(blk.mdd$adj.P.Val_jf22))
# consistent with the differences in the number of DEGs

# compare to us
resList <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results) {
	wt.mbv = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results, "_rev-gene-input_moderated-t-test.csv"))
	## Jaffe 2022
	### females first
	mbv.jf22_f = inner_join(filter(wt.mbv, coef=="F_NTC.MDD"),
		jf22[,c("gene_name","gene_id","logFC", "t", "P.Value", "adj.P.Val")],
		by="gene_id", suffix=c("_mbv","_bulk")) %>%
		mutate("MBv_coef"= "F_NTC.MDD", "MBv_annot"=results, bulk_source="Jaffe 2022") 

	tmp = filter(mbv.jf22_f, gene_name_mbv %in% ftest$V1)
	corList = list("corr_all" = list("jf22_f"=cor(mbv.jf22_f$logFC_mbv, mbv.jf22_f$logFC_bulk)),
		"corr_Ftest" = list("jf22_f"= cor(tmp$logFC_mbv, tmp$logFC_bulk))
		)

	### males next
	mbv.jf22_m = inner_join(filter(wt.mbv, coef=="M_NTC.MDD"),
	        jf22[,c("gene_name","gene_id","logFC", "t", "P.Value", "adj.P.Val")],
	        by="gene_id", suffix=c("_mbv","_bulk")) %>%
		mutate("MBv_coef"= "M_NTC.MDD",	"MBv_annot"=results, bulk_source="Jaffe 2022")
        corList[["corr_all"]][["jf22_m"]] <- cor(mbv.jf22_m$logFC_mbv, mbv.jf22_m$logFC_bulk)
                                          
        tmp = filter(mbv.jf22_m, gene_name_mbv %in% ftest$V1)
        corList[["corr_Ftest"]][["jf22_m"]] <- cor(tmp$logFC_mbv, tmp$logFC_bulk)

	## Girgenti 2021
	### females first
	mbv.gg21_f = inner_join(filter(wt.mbv, coef=="F_NTC.MDD"),
		gg21[,c("gene_id","gene_name","logFC","adj.P.Val")],
		by="gene_id", suffix=c("_mbv", "_bulk")) %>%
		mutate("MBv_coef"= "F_NTC.MDD",	"MBv_annot"=results, bulk_source="Girgenti 2021")

	corList[["corr_all"]][["gg21_f"]] <- cor(mbv.gg21_f$logFC_mbv, mbv.gg21_f$logFC_bulk)

	tmp = filter(mbv.gg21_f, gene_name_mbv %in% ftest$V1)
	corList[["corr_Ftest"]][["gg21_f"]] <- cor(tmp$logFC_mbv, tmp$logFC_bulk)

	### males next
	mbv.gg21_m = inner_join(filter(wt.mbv, coef=="M_NTC.MDD"),
	        gg21[,c("gene_id","gene_name","logFC","adj.P.Val")],
	        by="gene_id", suffix=c("_mbv", "_bulk")) %>%
		mutate("MBv_coef"= "M_NTC.MDD", "MBv_annot"=results, bulk_source="Girgenti 2021")

	corList[["corr_all"]][["gg21_m"]] <- cor(mbv.gg21_m$logFC_mbv, mbv.gg21_m$logFC_bulk)

	tmp = filter(mbv.gg21_m, gene_name_mbv %in% ftest$V1)
	corList[["corr_Ftest"]][["gg21_m"]] <- cor(tmp$logFC_mbv, tmp$logFC_bulk)

	# putting it all together
	tmp = bind_rows(mbv.jf22_f[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")],
		mbv.jf22_m[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")],
		mbv.gg21_f[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")],
		mbv.gg21_m[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")])

	tmp$MBv_coef = factor(tmp$MBv_coef, levels=c("F_NTC.MDD","M_NTC.MDD"))
	tmp$bulk_source = factor(tmp$bulk_source, levels=c("Jaffe 2022","Girgenti 2021"))

	return(list("data"=tmp, "logFC_correlations"=corList))
})

names(resList) <- c("domain.sp","domain.ct")

plot.df = bind_rows(filter(resList[[1]]$data, gene_name_mbv %in% ftest$V1),
		filter(resList[[2]]$data, gene_name_mbv %in% ftest$V1)) %>%
	mutate(MBv_annot= factor(MBv_annot, levels=c("smoothed-k9-1663","seurat-pc30"), labels=c("domain-SP","domain-CT")))

cor.df = data.frame(group= rep(names(resList[[1]]$logFC_correlations$corr_Ftest), 2),
	MBv_annot = c(rep("domain-SP", 4), rep("domain-CT", 4)),
	rho = c(unlist(resList[[1]]$logFC_correlations$corr_Ftest), unlist(resList[[2]]$logFC_correlations$corr_Ftest))
)
cor.df$MBv_coef = factor(cor.df$group, levels=names(resList[[1]]$logFC_correlations$corr_Ftest), labels=c("F_NTC.MDD","M_NTC.MDD","F_NTC.MDD","M_NTC.MDD"))
cor.df$bulk_source = factor(cor.df$group, levels=names(resList[[1]]$logFC_correlations$corr_Ftest), labels=c("Jaffe 2022", "Jaffe 2022", "Girgenti 2021", "Girgenti 2021"))

p1 <- ggplot(plot.df, aes(x=logFC_mbv, y=logFC_bulk))+
	geom_point(size=.3)+
	geom_smooth(method="lm", se=F, color="red3")+
	geom_text(data=cor.df, aes(x=-2, y=1, label=paste0("rho: ", round(rho, 3))), hjust=0, size=3, color="red3")+
	facet_grid(rows=vars(bulk_source), cols=vars(paste(MBv_coef, MBv_annot)), scales="free")+
#	coord_cartesian(ylim=c(-3,3), xlim=c(-3,3))+
	ggtitle("Whole-tissue model DE compared to bulk MDD")+
	theme_bw()+theme(panel.grid.minor=element_blank(), text=element_text(size=10))

# single BD option
xia24 = read.csv("raw-data/bulk_DE_results/Xia-2024_ST1-BD-DE-results.csv")
colnames(xia24)[1:2] = c("gene_id","gene_name")
colnames(xia24)[16:17] = c("logFC","adj.P.Val")

resList2 <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results) {
        wt.mbv = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results, "_rev-gene-input_moderated-t-test.csv"))

        ## Xia 2024
        ### females first
        mbv_f = inner_join(filter(wt.mbv, coef=="F_NTC.BPD"),
                xia24[,c("gene_id","gene_name","logFC","adj.P.Val")],
                by="gene_id", suffix=c("_mbv", "_bulk")) %>%
                mutate("MBv_coef"= "F_NTC.BPD", "MBv_annot"=results, bulk_source="Xia 2024")

        tmp = filter(mbv_f, gene_name_mbv %in% ftest$V1)
        corList = list("corr_all" = list("xia24_f"=cor(mbv_f$logFC_mbv, mbv_f$logFC_bulk)),
	        "corr_Ftest" = list("xia24_f"= cor(tmp$logFC_mbv, tmp$logFC_bulk))
                )

        ### males next
        mbv_m = inner_join(filter(wt.mbv, coef=="M_NTC.BPD"),
                xia24[,c("gene_id","gene_name","logFC","adj.P.Val")],
                by="gene_id", suffix=c("_mbv", "_bulk")) %>%
                mutate("MBv_coef"= "M_NTC.BPD", "MBv_annot"=results, bulk_source="Xia 2024")

        corList[["corr_all"]][["xia24_m"]] <- cor(mbv_m$logFC_mbv, mbv_m$logFC_bulk)

        tmp = filter(mbv_m, gene_name_mbv %in% ftest$V1)
        corList[["corr_Ftest"]][["xia24_m"]] <- cor(tmp$logFC_mbv, tmp$logFC_bulk)

        # putting it all together
        tmp = bind_rows(mbv_f[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")],
                mbv_m[,c("gene_id","gene_name_mbv","logFC_mbv","adj.P.Val_mbv","logFC_bulk","adj.P.Val_bulk","MBv_coef","MBv_annot","bulk_source")]
	)
        tmp$MBv_coef = factor(tmp$MBv_coef, levels=c("F_NTC.BPD","M_NTC.BPD"))

        return(list("data"=tmp, "logFC_correlations"=corList))
})
names(resList2) <- c("domain.sp","domain.ct")

plot.df2 = bind_rows(filter(resList2[[1]]$data, gene_name_mbv %in% ftest$V1),
                filter(resList2[[2]]$data, gene_name_mbv %in% ftest$V1)) %>%
        mutate(MBv_annot= factor(MBv_annot, levels=c("smoothed-k9-1663","seurat-pc30"), labels=c("domain-SP","domain-CT")))

cor.df2 = data.frame(group= rep(names(resList2[[1]]$logFC_correlations$corr_Ftest), 2),
        MBv_annot = c(rep("domain-SP", 2), rep("domain-CT", 2)),
        rho = c(unlist(resList2[[1]]$logFC_correlations$corr_Ftest), unlist(resList2[[2]]$logFC_correlations$corr_Ftest))
)
cor.df2$MBv_coef = factor(c("F_NTC.BPD","M_NTC.BPD","F_NTC.BPD","M_NTC.BPD"), levels=c("F_NTC.BPD","M_NTC.BPD"))
cor.df2$bulk_source = "Xia 2024"

p2 <- ggplot(plot.df2, aes(x=logFC_mbv, y=logFC_bulk))+
        geom_point(size=.3)+
        geom_smooth(method="lm", se=F, color="red3")+
        geom_text(data=cor.df2, aes(x=-2, y=1, label=paste0("rho: ", round(rho, 3))), hjust=0, size=3, color="red3")+
        facet_grid(rows=vars(bulk_source), cols=vars(paste(MBv_coef, MBv_annot)), scales="free")+
        #coord_cartesian(ylim=c(-3,3), xlim=c(-3,3))+
        ggtitle("Whole-tissue model DE compared to bulk BD")+
        theme_bw()+theme(panel.grid.minor=element_blank(), text=element_text(size=10))

ggsave(file="plots/08_compare_DEGs/bulk-RNAseq_whole-tissue-model.png",
	gridExtra::grid.arrange(p1, p2, layout_matrix=matrix(c(1,1,1,2,2), ncol=1)),
	width=7, height=6, bg="white")

