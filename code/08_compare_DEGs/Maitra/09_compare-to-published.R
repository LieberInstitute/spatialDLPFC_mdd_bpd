# exploratory comparison, run in interactive job and output not saved

library(dplyr)
library(ggplot2)

# load their results from published ST5
mt23 = read.csv("raw-data/Maitra-2023/Maitra-2023_all-DE-results.csv")

## is sig genes only, based on adj p value from within cluster-sex only (local)
mt23.sig = unique(mt23$gene)
length(mt23.sig) # 611 genes

## is union output from broad and cluster-specific
unique(mt23$cluster_id)


# load in my modeling of their data
mdd = read.csv("processed-data/08_compare_DEGs/Maitra/layer-restricted-pc2-ncells_low-res_moderated-t-test.csv")
mdd1 = filter(mdd, gene_name %in% mt23.sig)
length(unique(mdd1$gene_name)) # 541 genes they found to be sig were tested in our model


# start with astrocytes because clear and direct comparison possible given their use of a broad Ast category
# this comparison is a good place to start b/c they found only sig results in males, so its a simpler example to begin with
filter(mt23, cluster_id=="Ast") %>% group_by(Sex) %>% tally()
## 90 genes, all males

## first off, here is their results as a volcano, though missing the NS so hard to tell
ggplot(filter(mt23, cluster_id=="Ast"), aes(x=logFC, y=-log10(p_adj.loc)))+
	geom_point()
min(abs(filter(mt23, cluster_id=="Ast")$logFC)) # 0.83
## reasonable min abs logFC for sig results

mdd.ast = inner_join(filter(mdd1, cluster=="Ast") %>% select(logCPM=AveExpr, logFC, gene=gene_name, Sex=sex, p_val=P.Value, adj.P.Val), 
	filter(mt23, cluster_id=="Ast"), by=c("gene","Sex"), suffix=c("_mine","_theirs"))
table(mdd.ast$Sex)
## 84 of 90 genes were tested in our models

# checkout consistency in avg expr per group: it is meh
ggplot(mdd.ast, aes(x=logCPM_mine, y=logCPM_theirs))+
	geom_point()+
	# black line is a 1:1 match
	geom_abline(aes(slope=1, intercept=0))+
	# blue line is y~x fit
	geom_smooth(method="lm", se=F)+
	coord_cartesian(ylim=c(0,12), xlim=c(0,10))

p1 <- ggplot(mdd.ast, aes(x=logFC_mine, y=logFC_theirs))+
        geom_point(aes(color=-log10(p_val_mine)))+
        ggrepel::geom_text_repel(data=filter(mdd.ast, adj.P.Val<.05), aes(label=gene), min.segment.length=0, max.overlaps=Inf)+
	geom_hline(aes(yintercept=0), lty=3)+
	geom_vline(aes(xintercept=0), lty=3)+
	scale_color_viridis_c("-log10 raw\np-val (mine)", option="turbo")+
	labs(title="Ast - all M")

# now try inhib because also 1:1
filter(mt23, cluster_id=="InN") %>% group_by(Sex) %>% tally()  
## nevermind, only 1 gene

group_by(mt23, cluster_id) %>% tally() %>% arrange(desc(n))
#their highest sig groups are mostly cluster levels but also: Mic, OPC, L2/3
# but for L2/3 they have ExcN2 L2/3 and ExN9 L2/3

# Try Mic
filter(mt23, cluster_id=="Mic") %>% group_by(Sex) %>% tally()
## 74 genes, all F

ggplot(filter(mt23, cluster_id=="Mic"),	aes(x=logFC, y=-log10(p_adj.loc)))+
        geom_point()
min(abs(filter(mt23, cluster_id=="Mic")$logFC))	# 1.067

mdd.mic = inner_join(filter(mdd1, cluster=="M.V") %>% select(logCPM=AveExpr, logFC, gene=gene_name, Sex=sex, p_val=P.Value, adj.P.Val),
        filter(mt23, cluster_id=="Mic"), by=c("gene","Sex"), suffix=c("_mine","_theirs"))
table(mdd.mic$Sex)
## 62 of 74 genes were tested in our models

# waaaaay better consistency of avg logcounts
ggplot(mdd.mic, aes(x=logCPM_mine, y=logCPM_theirs))+
        geom_point()+
        # black line is a 1:1 match
        geom_abline(aes(slope=1, intercept=0))+	
        # blue line is y~x fit
        geom_smooth(method="lm", se=F)+
        coord_cartesian(ylim=c(0,12), xlim=c(0,12))

p2 <- ggplot(mdd.mic, aes(x=logFC_mine, y=logFC_theirs))+
        geom_point(aes(color=-log10(p_val_mine)))+
        ggrepel::geom_text_repel(data=filter(mdd.mic, adj.P.Val<.05), aes(label=gene), min.segment.length=0, max.overlaps=Inf)+
        geom_hline(aes(yintercept=0), lty=3)+
        geom_vline(aes(xintercept=0), lty=3)+
        scale_color_viridis_c("-log10 raw\np-val (mine)", option="turbo")+
        labs(title="Mic - all F")
# better than astro (not great, but better)

# try OPC
filter(mt23, cluster_id=="OPC") %>% group_by(Sex) %>% tally()
## 54 M, 3 F

ggplot(filter(mt23, cluster_id=="OPC"), aes(x=logFC, y=-log10(p_adj.loc)))+
        geom_point()
min(abs(filter(mt23, cluster_id=="OPC")$logFC)) # 0.535

mdd.opc = inner_join(filter(mdd1, cluster=="Olg") %>% select(logCPM=AveExpr, logFC, gene=gene_name, Sex=sex, p_val=P.Value, adj.P.Val),
        filter(mt23, cluster_id=="OPC"), by=c("gene","Sex"), suffix=c("_mine","_theirs"))
table(mdd.opc$Sex)
# 2 F, 51 M

# meh
ggplot(mdd.opc, aes(x=logCPM_mine, y=logCPM_theirs))+
        geom_point(aes(color=Sex))+
        # black line is a 1:1 match
        geom_abline(aes(slope=1, intercept=0))+ 
        # blue line is y~x fit
        geom_smooth(method="lm", se=F)+
        coord_cartesian(ylim=c(0,12), xlim=c(0,12))

p3 <- ggplot(mdd.opc, aes(x=logFC_mine, y=logFC_theirs))+
        geom_point(aes(color=-log10(p_val_mine)))+
        ggrepel::geom_text_repel(data=filter(mdd.opc, adj.P.Val<.05), aes(label=gene), min.segment.length=0, max.overlaps=Inf)+
        geom_hline(aes(yintercept=0), lty=3)+
        geom_vline(aes(xintercept=0), lty=3)+
        scale_color_viridis_c("-log10 raw\np-val (mine)", option="turbo")+
        labs(title="OPC - 51 M, 2 F")
# trend is better still tho 1 outlier


# try L2/3 ExN9
filter(mt23, cluster_id=="ExN9_L23") %>% group_by(Sex) %>% tally()
# 25 M, 6 F
ggplot(filter(mt23, cluster_id=="ExN9_L23"), aes(x=logFC, y=-log10(p_adj.loc)))+
        geom_point()
min(abs(filter(mt23, cluster_id=="ExN9_L23")$logFC)) # 0.69

mdd.exn9 = inner_join(filter(mdd1, cluster %in% c("L2","L3")) %>% select(logCPM=AveExpr, logFC, gene=gene_name, Sex=sex, p_val=P.Value, adj.P.Val, cluster_id=cluster),
        filter(mt23, cluster_id=="ExN9_L23"), by=c("gene","Sex"), suffix=c("_mine","_theirs"))
group_by(mdd.exn9, cluster_id_mine, Sex) %>% tally()
# 21 M, 4 F

ggplot(mdd.exn9, aes(x=logCPM_mine, y=logCPM_theirs))+
        geom_point(aes(color=Sex))+
	facet_grid(cols=vars(cluster_id_mine))+
        # black line is a 1:1 match
        geom_abline(aes(slope=1, intercept=0))+
        # blue line is y~x fit
        geom_smooth(method="lm", se=F)+
        coord_cartesian(ylim=c(0,12), xlim=c(0,12))

###### ALERT ######
###### this reminded me that the AveExpr is not per cluster but is overall in the fitted model for my results (identical L2 and L3)

p4 <- ggplot(mdd.exn9, aes(x=logFC_mine, y=logFC_theirs))+
        geom_point(aes(color=-log10(p_val_mine)))+facet_grid(cols=vars(cluster_id_mine))+
	ggrepel::geom_text_repel(data=filter(mdd.exn9, adj.P.Val<.05), aes(label=gene), min.segment.length=0, max.overlaps=Inf)+
        geom_hline(aes(yintercept=0), lty=3)+
        geom_vline(aes(xintercept=0), lty=3)+
        scale_color_viridis_c("-log10 raw\np-val (mine)", option="turbo")+
        labs(title="ExN9_L23 - 21 M, 4 F")

gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix=cbind(c(1,3,4),c(2,NA,4)))

#overall, not super consistent but also very different methods (and I didn't have age)
