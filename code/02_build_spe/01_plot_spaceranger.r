setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	#library(dplyr)
	library(ggplot2)
})
source("code/02_build_spe/getMBvSampleInfo_function.r")

l1 = list.files("processed-data/01_spaceranger")
summExists = file.exists(paste0("processed-data/01_spaceranger/",l1,"/outs/metrics_summary.csv"))
cat("\nFiles in 'processed-data/01_spaceranger' without summary csv:\n")
l1[summExists==FALSE]
l1 = l1[summExists]
cat("\nNumber of files in 'processed-data/01_spaceranger' with summary csv:\n")
length(l1)
sr = do.call(rbind, lapply(l1, function(x) read.csv(paste0("processed-data/01_spaceranger/",x,"/outs/metrics_summary.csv"))))

demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")
mdata = merge(sr, demo, by.x="Sample.ID", by.y="sample_id")

error.list = c(min(mdata$Valid.Barcodes)<.85, #changed from .9
	min(mdata$Sequencing.Saturation)<.8,
	min(mdata$Reads.Mapped.Confidently.to.Genome)<.7, #changed from .8
	max(mdata$Median.UMI.Counts.per.Spot)>8000) #changed from 5k

names(error.list) =c(paste("Min. value of Valid.Barcodes is",round(min(mdata$Valid.Barcodes),2),"-- change lower y limit"),
	paste("Min. value of Sequencing.Saturation is",round(min(mdata$Sequencing.Saturation),2),"-- change lower y limit"),
	paste("Min. value of Reads.Mapped.Confidently.to.Genome is",round(min(mdata$Reads.Mapped.Confidently.to.Genome),2),"-- change lower y limit"),
	paste("Max. value of Median.UMI.Counts.per.Spot is",round(max(mdata$Median.UMI.Counts.per.Spot),2),"-- change upper y limit"))

check.limits = sum(error.list)
if(sum(error.list)>0) {
	for(i in seq_along(error.list)) {
		if(isTRUE(error.list[[i]])) {
			warning(names(error.list)[i])}
	}
	stop("One or more variable requires adjusting y limits.")
}

#if want to plot everything from all three rounds separated by round
mdata$slide2 = paste(mdata$round, mdata$slide)
mdata$omit = ifelse(mdata$slide %in% c("V13Y10-020","V13B23-331"), TRUE, FALSE)

p1 <- ggplot(mdata, aes(x=slide2, y=Number.of.Reads))+
	geom_boxplot(aes(fill=omit))+scale_fill_manual(values=c("grey50","white"))+
	geom_text(aes(color=condition, label=brnum))+#default scale = ok
	labs(title="Number.of.Reads",x="seq. round & slide")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p2 <- ggplot(mdata, aes(x=slide2, y=Valid.Barcodes))+
	geom_boxplot(aes(fill=omit))+scale_fill_manual(values=c("grey50","white"))+
	geom_text(aes(color=condition, label=brnum))+ylim(.85,1)+
	labs(title="Valid.Barcodes", x="seq. round & slide")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p3 <- ggplot(mdata, aes(x=slide2, y=Sequencing.Saturation))+
	geom_boxplot(aes(fill=omit))+scale_fill_manual(values=c("grey50","white"))+
	geom_text(aes(color=condition, label=brnum))+ylim(.8,1)+
	labs(title="Sequencing.Saturation", x="seq. round & slide")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p4 <- ggplot(mdata, aes(x=slide2, y=Reads.Mapped.Confidently.to.Genome))+
	geom_boxplot(aes(fill=omit))+scale_fill_manual(values=c("grey50","white"))+
	geom_text(aes(color=condition, label=brnum))+ylim(.7,1)+
	labs(title="Reads.Mapped.Confidently.to.Genome", x="seq. round & slide")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p5 <- ggplot(mdata, aes(x=slide2, y=Median.UMI.Counts.per.Spot))+
	geom_boxplot(aes(fill=omit))+scale_fill_manual(values=c("grey50","white"))+
	geom_text(aes(color=condition, label=brnum))+ylim(0,8000)+
	labs(title="Median.UMI.Counts.per.Spot", x="seq. round & slide")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))

pdf("plots/02_build_spe/r1_r2_r3_spaceranger_overview.pdf", #height=12, width=12)
	height=6, width=36)
	#gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=2)
	p1
	p2
	p3
	p4
	p5
dev.off()
