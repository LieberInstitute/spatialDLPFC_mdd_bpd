setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(here)
})

l1 = list.files(here("processed-data","01_spaceranger"))
l1 = l1[file.exists(here("processed-data","01_spaceranger",l1,"outs","metrics_summary.csv"))]
sr = do.call(rbind, lapply(l1, function(x) read.csv(here("processed-data","01_spaceranger",x,"outs","metrics_summary.csv"))))
#demo = read.csv(here("raw-data","sample_info","MBv_demographics_combined.csv")) %>% mutate(sample_id=paste(slide, array, sep="_"))
#fix mbv demo table
demo = read.csv("raw-data/sample_info/MBv_demographics_combined.csv")
### first change 8 to B
demo$slide = gsub("V138","V13B",demo$slide)
### then change the 1 sample that had a different image number
demo[demo$mbv_sample=="MBv_081","slide"] = "V13B23-283"
demo[demo$mbv_sample=="MBv_081","array"] = "B1"
### then assign the samples missing slide info to the remaining slide
proc_slide = sapply(strsplit(l1,"_"), function(x) x[1])
missing_slide = setdiff(proc_slide, demo$slide)
demo[is.na(demo$slide),"slide"] = missing_slide

mdata = left_join(sr, mutate(demo, sample_id=paste(slide, array, sep="_")), by=c("Sample.ID"="sample_id"))

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

p1 <- ggplot(mdata, aes(x=slide, y=Number.of.Reads))+
	geom_boxplot()+
	geom_text(aes(color=condition, label=brain))+#default scale = ok
	labs(title="Number.of.Reads")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p2 <- ggplot(mdata, aes(x=slide, y=Valid.Barcodes))+
	geom_boxplot()+
	geom_text(aes(color=condition, label=brain))+ylim(.85,1)+
	labs(title="Valid.Barcodes")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p3 <- ggplot(mdata, aes(x=slide, y=Sequencing.Saturation))+
	geom_boxplot()+
	geom_text(aes(color=condition, label=brain))+ylim(.8,1)+
	labs(title="Sequencing.Saturation")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p4 <- ggplot(mdata, aes(x=slide, y=Reads.Mapped.Confidently.to.Genome))+
	geom_boxplot()+
	geom_text(aes(color=condition, label=brain))+ylim(.7,1)+
	labs(title="Reads.Mapped.Confidently.to.Genome")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
p5 <- ggplot(mdata, aes(x=slide, y=Median.UMI.Counts.per.Spot))+
	geom_boxplot()+
	geom_text(aes(color=condition, label=brain))+ylim(0,8000)+
	labs(title="Median.UMI.Counts.per.Spot")+
	theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1))
pdf(here("plots","02_build_spe","spaceranger_overview.pdf"), #height=12, width=12)
	height=6, width=36)
	#gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol=2)
	p1
	p2
	p3
	p4
	p5
dev.off()
