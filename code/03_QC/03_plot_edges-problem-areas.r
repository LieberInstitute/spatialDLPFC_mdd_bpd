setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(ggspavis)
	library(dplyr)
})

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
dim(spe)

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
#all edges are appropriately identified except for from the following samples
### V13B23-302_A1 and V13B23-302_C1
cdata$true_edges = ifelse(cdata$slide=="V13B23-302", FALSE, cdata$edge_outlier_genes)

tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>5, !is.na(problem_areas_genes.id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$problem_areas_genes.id)
length(unique(remove.areas))

stopifnot(identical(rownames(cdata),rownames(colData(spe))))

spe$edge_outlier_genes = cdata$edge_outlier_genes
spe$true_edges = cdata$true_edges
spe$problem_areas_genes.id = cdata$problem_areas_genes.id
spe$problem_areas_genes.size = cdata$problem_areas_genes.size

#problem area coldata
spe$edges_problem.areas_binary = ifelse(!is.na(spe$problem_areas_genes.id), "problem area", "none")
colData(spe)[spe$edge_outlier_genes, "edges_problem.areas_binary"] = "edge"

spe$problem_areas_grouped = ifelse(spe$problem_areas_genes.size<=5, "small", "flag")
spe$problem_areas_grouped = ifelse(spe$problem_areas_genes.size==0, "none", spe$problem_areas_grouped)
colData(spe)[spe$problem_areas_genes.id %in% remove.areas,"problem_areas_grouped"] = "remove"
colData(spe)[spe$true_edges,"problem_areas_grouped"] = "edge"

#plotting
spe2 <- spe[,spe$in_tissue]
dim(spe2)
spe2$edges_problem.areas_binary = as.factor(spe2$edges_problem.areas_binary)
table(spe2$edges_problem.areas_binary)
spe2$problem_areas_grouped = as.factor(spe2$problem_areas_grouped)
table(spe2$problem_areas_grouped)

spe2$dummy_slide = spe2$slide
colData(spe2)[spe2$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe2)[spe2$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

colData(spe2)$facet_violin = paste(spe2$round, spe2$dummy_slide)
colData(spe2)$facet_violin = ifelse(spe2$dummy_slide=="joint-283-339", "joint-283-339", colData(spe2)$facet_violin)
seed = levels(as.factor(colData(spe2)$facet_violin))

colData(spe2)$facet_spots = paste(spe2$round, spe2$sample_id)
colData(spe2)$facet_spots = ifelse(spe2$brnum=="Br5366", paste(spe2$facet_spots, spe2$brnum), spe2$facet_spots)

slideList2 = list(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList2 = lapply(slideList2, function(x) {
	unlist(lapply(x, function(y)
		sort(unique(colData(spe2)[spe2$facet_violin==y,"facet_spots"]))
  ))
})

color.palette = c("navy","#00a000","grey80")
names(color.palette) <- c("edge","problem area","none")

cat("\nGenerating raw results spot plots...\n")
plotList = lapply(slideList2, function(x) {
        l1 = x; names(l1) = x
        l1 = lapply(l1, function(y) spe2[,colData(spe2)$facet_spots==y])

        lapply(1:length(l1), function(z)
                suppressMessages(plotSpots(l1[[z]], annotate="edges_problem.areas_binary", #in_tissue=NULL,
                        point_size=0.2,
                        pal=color.palette)+
                        geom_point(show.legend=TRUE, size=.1)+
                        scale_color_manual("",values=color.palette, drop=F)+
                        labs(title=names(l1)[[z]])
                )
        )
})

pdf(file="plots/03_QC/edges-problem-areas_raw_spot-plots.pdf", width=12, height=16)
PRECAST::drawFigs(plotList[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
cat("\nPlot saved to: plots/03_QC/edges-problem-areas_raw_spot-plots.pdf\n")

#grouped results
color.palette2 = c("skyblue","palegoldenrod","#00a000","black","grey80")
names(color.palette2) <- c("edge","remove","flag","small","none")

cat("\nGenerating grouped spot plots...\n")
plotList2 = lapply(slideList2, function(x) {
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe2[,colData(spe2)$facet_spots==y])

	lapply(1:length(l1), function(z) 
		suppressMessages(plotSpots(l1[[z]], annotate="problem_areas_grouped", #in_tissue=NULL, 
			point_size=0.2,
			pal=color.palette2)+
			geom_point(show.legend=TRUE, size=.1)+
			scale_color_manual("",values=color.palette2, drop=F)+
			labs(title=names(l1)[[z]])
		)
	)
})


pdf(file="plots/03_QC/edges-problem-areas_grouped_spot-plots.pdf", width=12, height=16)
PRECAST::drawFigs(plotList2[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList2[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList2[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList2[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList2[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
cat("\nPlot saved to: plots/03_QC/edges-problem-areas_grouped_spot-plots.pdf\n")



cat("\nGenerating # low UMI bar plot...\n")
x_ordered = sort(unique(spe2$facet_spots))
colData(spe2)$x_facets = ""
colData(spe2)[spe2$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
colData(spe2)[spe2$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
colData(spe2)[spe2$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
colData(spe2)[spe2$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

tmp = mutate(as.data.frame(colData(spe2)), lowumi = sum_umi<=100, 
	     problem_areas_grouped=factor(problem_areas_grouped, levels=c("remove","flag","edge","small","none"))) %>%
  group_by(x_facets, facet_spots, problem_areas_grouped) %>%
  summarise(n_lowumi=sum(lowumi))

p2 <- ggplot(tmp, aes(x=facet_spots, y=n_lowumi, fill=problem_areas_grouped))+
  geom_bar(stat="identity", position="stack", color="black", linewidth=.5)+
  scale_y_continuous(expand=expansion(add=c(0,20)))+
  scale_fill_manual(values=color.palette2)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="# spots", fill="problem\nareas", title="Spots with <100 total UMI reads")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
ggsave(filename="plots/03_QC/edges-problem-areas_lowumi-spots_barplot.png", p2, bg="white", units="in", height=12, width=9)
cat("\nPlot saved to: plots/03_QC/edges-problem-areas_lowumi-spots_barplot.png\n")

cat("\nGenerating QC metric violin plots...\n")
spe2$problem_areas_grouped2 = ifelse(spe2$sum_umi<=100, "low UMI", "keep/flag")
spe2$problem_areas_grouped2 = ifelse(spe2$problem_areas_grouped %in% c("edge","remove"), "remove", spe2$problem_areas_grouped2)
spe2$problem_areas_grouped2 = as.factor(spe2$problem_areas_grouped2)
table(spe2$problem_areas_grouped2)

df = filter(as.data.frame(colData(spe2)), in_tissue==TRUE)
color.palette3 = c("palegoldenrod","grey80","red3")
names(color.palette3) = c("remove","keep/flag","low UMI")

p3 <- ggplot(df, aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=sum_umi, fill=problem_areas_grouped2))+
  geom_violin(position="dodge", scale="width")+
  scale_fill_manual(values=color.palette3)+
  geom_text(data=filter(df, problem_areas_grouped2=="low UMI") %>% group_by(condition, sex) %>% tally(),
            aes(label=paste(n, "spots"), y=200, fill=NULL), size=3, color="red3")+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(x="condition", y="sum_umi (log10 scale)", fill="", title="Library size")+
  theme_bw()

p4 <- ggplot(df, aes(x=factor(condition, levels=c("NTC","MDD","BPD")), 
                     y=sum_gene, fill=problem_areas_grouped2))+
  geom_violin(position="dodge", scale="width")+
  scale_fill_manual(values=color.palette3)+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:4)), labels=c("1","10","100","1k","10k"))+
  labs(x="condition", y="sum_gene (log10 scale)", fill="", title="Detected genes")+
  theme_bw()

p5 <- ggplot(mutate(df, expr_chrM_ratio=ifelse(is.na(expr_chrM_ratio), 0, expr_chrM_ratio)), 
                    aes(x=factor(condition, levels=c("NTC","MDD","BPD")), 
                     y=expr_chrM_ratio, fill=problem_areas_grouped2))+
  geom_violin(position="dodge", scale="width")+
  scale_fill_manual(values=color.palette3)+
  facet_wrap(vars(sex), ncol=1, scales="free_x")+
  labs(x="condition", y="expr_chrM_ratio", fill="", title="Mitochondrial fraction")+
  theme_bw()

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

tmp = group_by(df, sample_id, condition, sex, problem_areas_grouped2) %>% tally()
tmp2 = group_by(tmp, condition, sex, problem_areas_grouped2) %>% 
  mutate(outs1 = is_outlier(n)) %>% filter(outs1==T)

p6 <- ggplot(tmp, aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=n, 
                 fill=problem_areas_grouped2))+
  geom_boxplot(position="dodge")+facet_wrap(vars(sex))+
  geom_point(data=filter(tmp2, problem_areas_grouped2=="keep/flag"), 
             position=position_nudge(x=-.25), color="grey80", size=2, show.legend=F)+
  geom_point(data=filter(tmp2, problem_areas_grouped2=="low UMI"), 
             color="red3", size=2, show.legend=F)+
  geom_point(data=filter(tmp2, problem_areas_grouped2=="remove"), 
             position=position_nudge(x=.25), color="palegoldenrod", size=2, show.legend=F)+
  ggbreak::scale_y_break(c(1100,3000), scales="free")+
  scale_fill_manual(values=color.palette3)+
  labs(x="condition", y="# of spots", fill="", title="Spots passed QC (per sample)")+
  theme_bw()

ggsave(filename="plots/03_QC/edges-problem-areas_filtered-qc-metric_violin-plots.png", gridExtra::grid.arrange(p3, p4, p5, print(p6), ncol=2), 
       bg="white", units="in", height=9, width=12)
cat("\nPlot saved to: plots/03_QC/edges-problem-areas_filtered-qc-metric_violin-plots.png\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
