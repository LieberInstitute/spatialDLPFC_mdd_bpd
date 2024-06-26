setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggplot2)
	library(dplyr)
	library(here)
})

l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[grep("^bindev_V.{9}_default-brain",l1)]

bindev.df = do.call(rbind, lapply(l1, function(x) {
		tmp = read.csv(here("processed-data","04_preprocessing",x))
		y=substr(x,8,17)
		tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
		return(tmp)
	})
)

#limit search for biased features to top 10k ranked genes
seq1 = seq(0,10000, by=1000)
bindev.df$rank_default_bin = NA
for (i in 2:length(seq1)) {
	bindev.df = mutate(bindev.df, 
		rank_default_bin=if_else(between(rank_default, seq1[i-1]+1, seq1[i]), 
			paste0("[",seq1[i-1]+1,",",seq1[i],"]"),
			rank_default_bin))
}
bindev.df = filter(bindev.df, !is.na(rank_default_bin))

#determine pooled SD
mean1 = mean(bindev.df$r.diff)
sd1 = sd(bindev.df$r.diff)
cat("\npooled mean rank diff.=",mean1)
cat("\npooled SD rank diff.=",sd1,"\n")
bindev.df$nSD = (bindev.df$r.diff-mean1)/sd1
bindev.df$nSD.bin = cut(abs(bindev.df$nSD), right=FALSE, breaks=seq(0,max(bindev.df$nSD)+5, by=5), include.lowest=TRUE)

#determine per slide SD
cat("\nper-slide mean and SD rank diff.\n")
group_by(bindev.df, slide) %>% summarise(avg.r.diff=mean(r.diff), sd.r.diff=sd(r.diff))
bindev.df = group_by(bindev.df, slide) %>% mutate(slide.mean=mean(r.diff), slide.sd=sd(r.diff), nSD.slide=(r.diff-slide.mean)/slide.sd) %>% ungroup()
bindev.df$nSD.slide.bin = cut(abs(bindev.df$nSD.slide), right=FALSE, breaks=seq(0,max(bindev.df$nSD.slide)+5, by=5), include.lowest=TRUE)

#id outliers at 5SD threshold
cat("\nsubject-biased genes identified at nSD(rank diff.)>=5\n")
bindev.df = mutate(ungroup(bindev.df), nSD.outlier=nSD.bin!="[0,5)", nSD.slide.outlier=nSD.slide.bin!="[0,5)",
	outlier.group=factor(paste(nSD.outlier, nSD.slide.outlier),
		levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
		labels=c("none","per slide only","pooled only","both")
		)
	)

#save dataframe
write.csv(bindev.df, here("processed-data","04_preprocessing","bindev_biased-feature_results_10k.csv"), row.names=F)
#save biased gene list
outlier.list = list("nSD5.pooled.outlier"=unique(filter(bindev.df, nSD.bin!="[0,5)")$gene),
	"nSD5.slide.outlier"=unique(filter(bindev.df, nSD.slide.bin!="[0,5)")$gene))
saveRDS(outlier.list, here("processed-data","04_preprocessing","bindev_biased-feature_5sd_list.rda"))

#plots!
#density plot faceted by slide showing that the distribution of rank differences is stable with increasing default rank (up to 10k)
col.pal = RColorBrewer::brewer.pal(n=10,"RdYlBu")
p1 <- ggplot(bindev.df, aes(r.diff, color=rank_default_bin))+
	geom_density()+
	facet_wrap(vars(slide), ncol=2)+
	scale_color_manual(values=col.pal[10:1])+
	scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1),
		breaks=c(-100,-10,-1,10^(0:3)), labels=format(c(-100,-10,-1,10^(0:3)), scientific=F))+
	labs(x="change in rank with batch consideration", title="batch = subject", fill="rank (no batch)")+
	theme_bw()

ggsave(file=here("plots","04_preprocessing","biased-features-10k_rank-diff_density.png"), plot=p1,
        height=8, width=6, bg="white")
cat("\nfaceted density plot save to:",here("plots","04_preprocessing","biased-features-10k_rank-diff_density.png"))
#histogram of nSD for both threshold methods
p2 <- ggplot(bindev.df, aes(x=abs(nSD), fill=nSD.bin))+
	geom_histogram(color="grey90")+
	scale_fill_brewer(palette="YlOrRd")+
	facet_wrap(vars(slide), ncol=2)+
	labs(title="pooled SD calculation", x="abs(nSD)",fill="nSD",y="# genes")+
	scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
		breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
	theme_bw()
p3 <- ggplot(bindev.df, aes(x=abs(nSD.slide), fill=nSD.slide.bin))+
	geom_histogram(color="grey90")+
	scale_fill_brewer(palette="YlOrRd")+
	facet_wrap(vars(slide), ncol=2)+
	labs(title="per-slide SD calculation", x="abs(nSD)",fill="nSD",y="# genes")+
	scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
		breaks=10^(0:4), labels=format(10^(0:4), scientific=F))+
	theme_bw()

#rank scatter plot highlighting differences between two approaches
p4 <- ggplot(bindev.df, aes(x=rank_default, y=rank_brain, color=outlier.group))+
	geom_point(size=.5)+facet_wrap(vars(slide), ncol=2)+
	labs(x="rank (no batch)",y="rank (batch= subject)", color="outlier",title="candidate subject-biased genes")+
	scale_color_manual(values=c("grey","red","blue","black"))+
	theme_bw()

ggsave(file=here("plots","04_preprocessing","biased-features-10k_pooled-sd_per-slide-sd_hist-and-scatter.png"), plot=gridExtra::grid.arrange(p2,p3,p4, ncol=3),
        height=8, width=20, bg="white")
cat("\nthree figure pooled and per-slide SD threshold plot saved to:",here("plots","04_preprocessing","biased-features-10k_pooled-sd_per-slide-sd_hist-and-scatter.png"))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
