# tested in interactive node, need to integrate with rest of workflow
library(dplyr)
library(ggplot2)

ftest = read.table("raw-data/SCENIC_aux/tf_lists/MBv_PRECAST-Seurat_F-test-adjp-05.txt")

ch23 = read.csv("raw-data/bulk_DE_results/Chatzinakos-2023_ST6B-MDD-DE.csv")
colnames(ch23)[6:7] = c("adj.P.Val","gene_name")

ch23.bulk = filter(ch23, Cell_type_Bulk=="Bulk")

resList <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results) {
	wt.mbv = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results, "_rev-gene-input_moderated-t-test.csv"))

	### females first
	mdd_f = inner_join(filter(wt.mbv, coef=="F_NTC.MDD"),
		ch23.bulk[,c("gene_name","logFC", "t", "P.Value", "adj.P.Val")],
		by="gene_name", suffix=c("_mbv","_bulk")) %>%
		mutate("MBv_coef"= "F_NTC.MDD", "MBv_annot"=results, bulk_source="Chatzinakos 2023") 

	tmp = filter(mdd_f, gene_name %in% ftest$V1)
	corList = list("corr_all" = list("ch23_f"=cor(mdd_f $logFC_mbv, mdd_f $logFC_bulk)),
		"corr_Ftest" = list("ch23_f"= cor(tmp$logFC_mbv, tmp$logFC_bulk))
		)

	### males next
	mdd_m = inner_join(filter(wt.mbv, coef=="M_NTC.MDD"),
	        ch23.bulk[,c("gene_name","logFC", "t", "P.Value", "adj.P.Val")],
	        by="gene_name", suffix=c("_mbv","_bulk")) %>%
		mutate("MBv_coef"= "M_NTC.MDD",	"MBv_annot"=results, bulk_source="Chatzinakos 2023")
        corList[["corr_all"]][["ch23_m"]] <- cor(mdd_m $logFC_mbv, mdd_m $logFC_bulk)
                                          
        tmp = filter(mdd_m, gene_name %in% ftest$V1)
        corList[["corr_Ftest"]][["ch23_m"]] <- cor(tmp$logFC_mbv, tmp$logFC_bulk)

        # putting it all together
        tmp = bind_rows(mdd_f, mdd_m)
        tmp$MBv_coef = factor(tmp$MBv_coef, levels=c("F_NTC.MDD","M_NTC.MDD"))

        return(list("data"=tmp, "logFC_correlations"=corList))
})
names(resList) <- c("domain.sp","domain.ct")

plot.df = bind_rows(filter(resList[[1]]$data, gene_name %in% ftest$V1),
                filter(resList[[2]]$data, gene_name %in% ftest$V1)) %>%
#	bind_rows(resList[[1]]$data, resList[[2]]$data) %>%
        mutate(MBv_annot= factor(MBv_annot, levels=c("smoothed-k9-1663","seurat-pc30"), labels=c("domain-SP","domain-CT")))

cor.df = data.frame(group= rep(names(resList[[1]]$logFC_correlations$corr_Ftest), 2),
        MBv_annot = c(rep("domain-SP", 2), rep("domain-CT", 2)),
        rho = c(unlist(resList[[1]]$logFC_correlations$corr_Ftest), unlist(resList[[2]]$logFC_correlations$corr_Ftest))
#        rho = c(unlist(resList[[1]]$logFC_correlations$corr_all), unlist(resList[[2]]$logFC_correlations$corr_all))
)
cor.df$MBv_coef = factor(c("F_NTC.MDD","M_NTC.MDD","F_NTC.MDD","M_NTC.MDD"), levels=c("F_NTC.MDD","M_NTC.MDD"))
cor.df$bulk_source = "Chatzinakos 2023"

p <- ggplot(plot.df, aes(x=logFC_mbv, y=logFC_bulk))+
        geom_point(size=.3)+
        geom_smooth(method="lm", se=F, color="red3")+
        geom_text(data=cor.df, aes(x=-2.5, y=4, label=paste0("rho: ", round(rho, 3))), hjust=0, size=3, color="red3")+
        facet_grid(rows=vars(bulk_source), cols=vars(paste(MBv_coef, MBv_annot)))+#, scales="free")+
        coord_cartesian(ylim=c(-5,5), xlim=c(-3,3))+
        labs(title="Whole-tissue model DE compared to similar sn MDD model", subtitle="(snRNA-seq is sex-pooled with covariate)",
		y="logFC_sn")+
        theme_bw()+theme(panel.grid.minor=element_blank(), text=element_text(size=10))
ggsave("plots/08_compare_DEGs/sn-Chatz_whole-tissue-model_F-sig-genes.png",
	p, width=7, height=3, bg="white")

# layer restricted comparison
mbv = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_moderated-t-test.csv")
mbv$cluster = factor(mbv$cluster, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"),
	labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))
mbv$source = "MBv"

## just domain-CT F-test for this
ftest = read.csv("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_rev-gene-input_F-test_all.csv") %>% filter(adj.P.Val<.05)
colnames(ftest)[56] = "V1"

# astro first
mbv.ast = filter(mbv, cluster=="Ast", group=="NTC.MDD")
ch23.ast = filter(ch23, Cell_type_Bulk=="Astro") %>% rename("cluster" = Cell_type_Bulk) %>% mutate(cluster="Ast")
ast = inner_join(mbv.ast[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster","sex","group")], 
	ch23.ast[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster")],
	by="gene_name", suffix=c("_mbv","_sn"))
max(abs(ast$t_mbv))
max(abs(ast$t_sn))

ggplot(ast, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	facet_wrap(vars(sex))+
	coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
	labs(x="sex-specific MBv t stat", y="sex-pooled Chatz. t stat")+
	theme_bw()


# inhib next
mbv.inb = filter(mbv, cluster=="Inb", group=="NTC.MDD")
ch23.inb = filter(ch23, Cell_type_Bulk=="In") %>% rename("cluster" = Cell_type_Bulk) %>% mutate(cluster="Inb")
inb = inner_join(mbv.inb[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster","sex","group")], 
	ch23.inb[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster")],
	by="gene_name", suffix=c("_mbv","_sn"))
max(abs(inb$t_mbv))
max(abs(inb$t_sn))

ggplot(inb, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	facet_wrap(vars(sex))+
	coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
	labs(x="sex-specific MBv t stat", y="sex-pooled Chatz. t stat")+
	theme_bw()

# olg
mbv.olg = filter(mbv, cluster=="Olg", group=="NTC.MDD")
mbv.olg = bind_rows(mutate(mbv.olg, cluster_match="Oligo"), mutate(mbv.olg, cluster_match="OPC"))
ch23.olg = filter(ch23, Cell_type_Bulk %in% c("Oligo","OPC")) %>% rename("cluster_match" = Cell_type_Bulk)
olg = inner_join(mbv.olg[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster","cluster_match","sex","group")], 
	ch23.olg[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster_match")],
	by=c("gene_name","cluster_match"), suffix=c("_mbv","_sn"))
max(abs(olg$t_mbv))
max(abs(olg$t_sn))

ggplot(olg, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	facet_grid(cols=vars(sex), rows=vars(cluster_match))+
	coord_cartesian(xlim=c(-8,8), ylim=c(-8,8))+
	labs(x="sex-specific MBv t stat", y="sex-pooled Chatz. t stat")+
	theme_bw()


# mv
mbv.mv = filter(mbv, cluster=="M.V", group=="NTC.MDD")
mbv.mv = bind_rows(mutate(mbv.mv, cluster_match="Endo"), mutate(mbv.olg, cluster_match="Per"), mutate(mbv.olg, cluster_match="Microglia"))
ch23.mv = filter(ch23, Cell_type_Bulk %in% c("Endo","Per","Microglia")) %>% rename("cluster_match" = Cell_type_Bulk)
mv = inner_join(mbv.mv[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster","cluster_match","sex","group")], 
	ch23.mv[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster_match")],
	by=c("gene_name","cluster_match"), suffix=c("_mbv","_sn"))
max(abs(mv$t_mbv))
max(abs(mv$t_sn))

ggplot(mv, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	facet_grid(cols=vars(sex), rows=vars(cluster_match))+
	coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
	labs(x="sex-specific MBv t stat", y="sex-pooled Chatz. t stat")+
	theme_bw()

most.df = bind_rows(mutate(ast, cluster_match="Ast"), mutate(inb, cluster_match="Inb"),
	olg, mv) %>%
	mutate(cluster_match = factor(cluster_match, levels=c("Ast","Inb","Oligo","OPC","Endo","Per","Microglia")))

cor.df = do.call(rbind, lapply(levels(most.df$cluster_match), function(x) {
	f1 = filter(most.df, cluster_match==x, sex=="F")
	m1 = filter(most.df, cluster_match==x, sex=="M")
	f2 = filter(f1, gene_name %in% ftest$V1)
	m2 = filter(m1, gene_name %in% ftest$V1)
	data.frame("cluster_match"=rep(x,4),
		sex=rep(c("F","M"),2),
		stat=rep("t",4),
		gene_set=c(rep("all",2),rep("F-sig",2)),
		rho=c(cor(f1$t_mbv, f1$t_sn), cor(m1$t_mbv, m1$t_sn),
			cor(f2$t_mbv, f2$t_sn), cor(m2$t_mbv, m2$t_sn))
	)
})) %>% mutate(cluster_match= factor(cluster_match, levels=levels(most.df$cluster_match)))


ggplot(most.df, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	facet_grid(cols=vars(sex), rows=vars(cluster_match))+
	coord_cartesian(xlim=c(-8,8), ylim=c(-8,8))+
	labs(x="sex-specific MBv t stat", y="sex-pooled Chatz. t stat")+
	theme_bw()

filt.df = filter(most.df, gene_name %in% ftest$V1)
max(abs(filt.df$t_mbv))
max(abs(filt.df$t_sn))
p1 <- ggplot(filt.df, aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
	geom_smooth(method="lm", color="red3", se=F)+
	geom_text(data=filter(cor.df, gene_set=="F-sig"), aes(x=-7, y=5, label=paste0("rho: ", round(rho, 3))), hjust=0, size=3, color="red3")+
	facet_grid(cols=vars(sex), rows=vars(cluster_match))+
	coord_cartesian(xlim=c(-7.5,7.5), ylim=c(-7.5,7.5))+
	labs(x="t stat (MBv SRT)", y="t stat (snRNA-seq)",
		title="Domain-CT restricted model DE vs pub",
		subtitle="Chatzinakos 2023: sex-pooled MDD (facets)\n1 domain-CT to 1 or multi cell type(s)")+
	theme_bw()+theme(plot.title.position="plot")
ggsave("plots/08_compare_DEGs/sn-Chatz_domain-restricted-model_part1_domain-ct-F-sig-genes.png",
	p1, height= 10, width=4, bg="white")

# ex
mbv.ex = filter(mbv, cluster %in% c("L2.3","L4","L5","L6"), group=="NTC.MDD") %>% mutate(cluster_match="Ex")
ch23.ex = filter(ch23, Cell_type_Bulk=="Ex") %>% rename("cluster_match" = Cell_type_Bulk)
ex = inner_join(mbv.ex[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster","cluster_match","sex","group")], 
	ch23.ex[,c("logFC","t","P.Value","adj.P.Val","gene_name","cluster_match")],
	by=c("gene_name","cluster_match"), suffix=c("_mbv","_sn"))
max(abs(ex$t_mbv))
max(abs(ex$t_sn))

ex$cluster = droplevels(ex$cluster)
cor.df2 = do.call(rbind, lapply(levels(ex$cluster), function(x) {
        f1 = filter(ex, cluster==x, sex=="F")
        m1 = filter(ex, cluster==x, sex=="M")
        f2 = filter(f1, gene_name %in% ftest$V1)
        m2 = filter(m1, gene_name %in% ftest$V1)
        data.frame("cluster"=rep(x,4),
                sex=rep(c("F","M"),2),
                stat=rep("t",4),
                gene_set=c(rep("all",2),rep("F-sig",2)),
                rho=c(cor(f1$t_mbv, f1$t_sn), cor(m1$t_mbv, m1$t_sn),
                        cor(f2$t_mbv, f2$t_sn), cor(m2$t_mbv, m2$t_sn))
        )
}))

p2 <- ggplot(filter(ex, gene_name %in% ftest$V1), aes(x=t_mbv, y=t_sn))+
	geom_point(size=.1)+
        geom_smooth(method="lm", color="red3", se=F)+
        geom_text(data=filter(cor.df2, gene_set=="F-sig"), aes(x=-7, y=5, label=paste0("rho: ", round(rho, 3))), hjust=0, size=3, color="red3")+
	facet_grid(cols=vars(cluster), rows=vars(sex))+
	coord_cartesian(xlim=c(-7,7), ylim=c(-7,7))+
	labs(x="t stat (MBv SRT)", y="t stat (snRNA-seq)",
		title="Domain-CT restricted model DE vs published results",
		subtitle="Chatzinakos 2023: sex-pooled MDD\nSingle ExcN cell type (y-axis) vs my ExcN domains-CT")+
	theme_bw()
ggsave("plots/08_compare_DEGs/sn-Chatz_domain-restricted-model_part2_domain-ct-F-sig-genes.png",
        p2, height= 4,	width=8, bg="white")
