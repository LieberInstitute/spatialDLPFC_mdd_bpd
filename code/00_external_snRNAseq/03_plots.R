library(dplyr)
library(ggplot2)

sn.bd = read.csv("processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_filtered_observations.csv") %>%
	filter(individualID!="BD24")
sn.mdd = read.csv("processed-data/00_external_snRNAseq/Maitra2023_control-MDD_filtered_observations.csv")

# detected genes
ggplot(sn.bd, aes(x=nFeature_RNA, color=individualID))+
	stat_ecdf()+coord_cartesian(xlim=c(0,10000))+
	facet_grid(cols=vars(Biological_Sex), rows=vars(Disorder))+
	theme(legend.position="none")

ggplot(sn.mdd, aes(x=nFeature_RNA, color=Sample))+
	stat_ecdf()+coord_cartesian(xlim=c(0,10000))+
	facet_grid(cols=vars(Sex), rows=vars(Condition))+
	theme(legend.position="none")

# library size
ggplot(sn.bd, aes(x=nCount_RNA, color=individualID))+
        stat_ecdf()+coord_cartesian(xlim=c(0,20000))+
        facet_grid(cols=vars(Biological_Sex), rows=vars(Disorder))+
        theme(legend.position="none")

ggplot(sn.mdd, aes(x=nCount_RNA, color=Sample))+
        stat_ecdf()+coord_cartesian(xlim=c(0,20000))+
        facet_grid(cols=vars(Sex), rows=vars(Condition))+
        theme(legend.position="none")

# versus
ggplot(sn.bd, aes(x=log10(nFeature_RNA), y=log10(nCount_RNA)))+
	geom_point(size=.1)+
	facet_grid(rows=vars(Disorder), cols=vars(Biological_Sex))

ggplot(sn.mdd, aes(x=log10(nFeature_RNA), y=log10(nCount_RNA)))+
        geom_point(size=.1)+
        facet_grid(rows=vars(Condition), cols=vars(Sex))

## sn.mdd by chemistry and sequencing tech
ggplot(sn.mdd, aes(x=log10(nFeature_RNA), y=log10(nCount_RNA), color=Sex))+
	geom_point(size=.1)+
	facet_grid(rows=vars(Chemistry), cols=vars(Sequencing))


# nuclei per donor
d.bd = group_by(sn.bd, individualID, Disorder, Biological_Sex) %>% tally() %>%
        mutate(cond_sex= factor(paste(factor(Disorder, levels=c("control","Bipolar Disorder"), labels=c("con","bd")),
                        factor(Biological_Sex, levels=c("female","male"), labels=c("F","M"))),
                levels=c("con F","con M","bd F","bd M")))
d.mdd = group_by(sn.mdd, Sample, Condition, Sex) %>% tally() %>%
	mutate(cond_sex= factor(paste(factor(Condition, levels=c("Control","Case"), labels=c("con","mdd")),
			factor(Sex, levels=c("Female","Male"), labels=c("F","M"))),
		levels=c("con F","con M","mdd F","mdd M")))
donor.df = bind_rows(mutate(d.bd[,c("n","cond_sex")], source="SZBD", outlier=n>18000),
	mutate(d.mdd[,c("n","cond_sex")], source="Maitra", outlier=F)) %>%
	mutate(source= factor(source, levels=c("SZBD","Maitra")),
		outlier= factor(as.character(outlier), levels=c("TRUE","FALSE"), labels=c("BD4","none")))

p1 <- ggplot(donor.df, aes(x=cond_sex, y=n, color=outlier))+
	ggbeeswarm::geom_beeswarm()+scale_color_manual(values=c("none"="black", "BD4"="red3"))+
	facet_grid(cols=vars(source), scales="free_x")+
	ylim(0,20000)+labs(title="Total nuclei per donor", y="# nuclei", x="dx-sex group")+
	theme_bw()

# spots per donor
mdata = read.csv("processed-data/93_globus/MBv_n119_observations.csv.gz", row.names=1)
d.mbv1 = filter(mdata, !is.na(domain_sp)) %>% group_by(brnum, diagnosis, sex) %>% tally() %>%
	mutate(cond_sex= factor(paste(diagnosis, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BD F","BD M")))
d.mbv2 = filter(mdata, !is.na(domain_ct)) %>% group_by(brnum, diagnosis, sex) %>% tally() %>%
        mutate(cond_sex= factor(paste(diagnosis, sex), levels=c("NTC F","NTC M","MDD F","MDD M","BD F","BD M")))
donor.df2 = bind_rows(mutate(d.mbv1, source="MBv (domain-SP)"),
	mutate(d.mbv2, source="MBv (domain-CT)")) %>%
	mutate(source= factor(source, levels=c("MBv (domain-SP)", "MBv (domain-CT)")))

p2 <- ggplot(donor.df2, aes(x=cond_sex, y=n))+
        ggbeeswarm::geom_beeswarm()+
        facet_grid(cols=vars(source), scales="free_x")+
        ylim(0,5500)+labs(title="Total spots per donor", y="# spots", x="dx-sex group")+
        theme_bw()

# donor per dx-sex group
group.df = bind_rows(group_by(donor.df, source, cond_sex, outlier) %>% tally(),
        group_by(donor.df2, source, cond_sex, outlier="none") %>% tally()) %>%
        mutate(outlier=	factor(as.character(outlier), levels=c("BD4","none")))

p3 <- ggplot(group.df, aes(x=cond_sex, y=n, fill=outlier))+
        geom_bar(stat="identity", position="stack")+scale_fill_manual(values=c("none"="black", "BD4"="red3"))+
        facet_grid(cols=vars(source), scales="free_x", space="free_x")+
        ylim(0,25)+
        labs(title="Donors per dx-sex group", y="# donors", x="dx-sex group")+
        theme_bw()

# observations per dx-sex group
group.df2 = bind_rows(group_by(donor.df, source, cond_sex, outlier) %>% summarise(n_group=sum(n)),
	group_by(donor.df2, source, cond_sex, outlier="none") %>% summarise(n_group=sum(n))) %>%
	mutate(outlier= factor(as.character(outlier), levels=c("BD4","none")))

p4 <- ggplot(group.df2, aes(x=cond_sex, y=n_group, fill=outlier))+
        geom_bar(stat="identity", position="stack")+scale_fill_manual(values=c("none"="black", "BD4"="red3"))+
        facet_grid(cols=vars(source), scales="free_x", space="free_x")+
        ylim(0,100000)+
	labs(title="Total observations per dx-sex group", y="# nuclei or spots", x="dx-sex group")+
        theme_bw()

gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix=rbind(c(1,2),c(3,3),c(4,4)))

# feature lists
gn.bd = read.csv("processed-data/00_external_snRNAseq/SZBDMultiseq_control-BD_unfiltered_features.csv")
gn.mdd = read.csv("processed-data/00_external_snRNAseq/Maitra2023_control-MDD_unfiltered_features.csv")

rdata = read.csv("processed-data/93_globus/MBv_n119_features.csv.gz", row.names=1)

tmp = inner_join(gn.mdd, gn.bd, by=c("gene_id"="featureid"))
tmp = inner_join(tmp, rdata[,c("gene_id","gene_name","n_spots")])
nrow(tmp) # 34291

## between snRNAseq correlation
cor(tmp$n_nuclei, tmp$n_cells) # 0.97542
## SRT to snRNAseq correlations
cor(tmp$n_spots, tmp$n_nuclei) # 0.658491
cor(tmp$n_spots, tmp$n_cells) # 0.6095497

# recheck correlation when filtered to genes that passed minimum obs filters
tmp2 = filter(tmp, n_nuclei>100, n_cells>50, n_spots>5)
nrow(tmp2) # 23801

## between snRNAseq correlation
cor(tmp2$n_nuclei, tmp2$n_cells) # 0.9719848 
## SRT to snRNAseq correlations
cor(tmp2$n_spots, tmp2$n_nuclei) # 0.6227164
cor(tmp2$n_spots, tmp2$n_cells) # 0.5658171

# check when filters are to DE genes
de.genes = rownames(rdata)[rdata$DE_domain.sp | rdata$DE_domain.ct]
length(de.genes) # 13948

length(setdiff(de.genes, gn.bd$featureid)) # 1010 genes in DE input that were removed from SZBD prior to synapse upload
length(setdiff(de.genes, gn.mdd$gene_id)) # 0

tmp3 = filter(tmp, gene_id %in% de.genes)
dim(tmp3) # 12938     7

## between snRNAseq correlation
cor(tmp3$n_nuclei, tmp3$n_cells) # 0.9623361
## SRT to snRNAseq correlations
cor(tmp3$n_spots, tmp3$n_nuclei) # 0.5059636
cor(tmp3$n_spots, tmp3$n_cells) # 0.4204839
