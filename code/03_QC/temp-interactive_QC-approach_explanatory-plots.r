suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(ggplot2)
  library(dplyr)
})

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
dim(spe)
tmp = spe[,spe$sample_id=="V13Y10-023_A1"]
tmp$lg10.umi = log10(tmp$sum_umi)
colSums(is.na(colData(tmp)))
tmp$lg10.genes = log10(tmp$sum_gene)
colSums(is.na(colData(tmp)))

ggspavis::plotSpots(tmp, annotate="lg10.umi", point_size=1)+
  scale_color_gradient(low="white", high="navy")+#, labels=function(x) paste0(x/1000,"k"))+
  labs(title="r1 V13Y10-023_A1", color="log10(UMI)")+
  theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
ggspavis::plotSpots(tmp, annotate="lg10.genes", point_size=1)+
  scale_color_gradient(low="white", high="navy")+#, labels=function(x) paste0(x/1000,"k"))+
  labs(title="r1 V13Y10-023_A1", color="log10(genes)")+
  theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
filter(cdata, in_tissue==TRUE) %>% group_by(sample_id) %>% summarise(n_spots_removed=sum(remove_spots))


length(unique(cdata$problem_areas_genes.id)) #4069
filter(cdata, problem_areas_genes.size>5) %>% distinct(problem_areas_genes.id) %>% nrow() #443
filter(cdata, problem_areas_genes.size>10) %>% distinct(problem_areas_genes.id) %>% nrow() #279
filter(cdata, problem_areas_genes.size>20) %>% distinct(problem_areas_genes.id) %>% nrow() #181
filter(cdata, problem_areas_genes.size>50) %>% distinct(problem_areas_genes.id) %>% nrow() #90

cdata$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size<=5, "<=5", "flag")
cdata$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size==0, "none", cdata$problem_areas_grouped)
cdata$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size>5, "(5,10]", cdata$problem_areas_grouped)
cdata$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size>10, "(10,20]", cdata$problem_areas_grouped)
cdata$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size>20, ">20", cdata$problem_areas_grouped)
cdata$problem_areas_grouped = ifelse(cdata$true_edges==TRUE, "edge", cdata$problem_areas_grouped)
cdata$problem_areas_grouped = factor(cdata$problem_areas_grouped, levels=c("edge",">20","(10,20]","(5,10]","<=5","none"))
table(cdata$problem_areas_grouped)
table(distinct(cdata, problem_areas_genes.id, problem_areas_grouped)$problem_areas_grouped)

#cp1 = c("black","grey50","grey","white")
cp1 = c("grey90","grey","grey50","black")
names(cp1) <- c(">20","(10,20]","(5,10]","<=5")
p1 <- ggplot(filter(cdata, in_tissue==TRUE, !problem_areas_grouped %in% c("none","edge")) %>% group_by(problem_areas_grouped, problem_areas_genes.id) %>% tally(),
       aes(x=n, fill=problem_areas_grouped))+
  geom_histogram(bins=60, color="black", linewidth=.3)+scale_x_log10()+scale_y_log10()+
  #ggbreak::scale_y_break(c(400,2500))+
  scale_fill_manual(values=cp1)+
  theme_bw()+
  labs(x="problem area size (# spots)", y="# of problem areas", 
       title="Non-edge problem areas", fill="problem\narea size")


p2 <- ggplot(filter(cdata, in_tissue==TRUE) %>% group_by(problem_areas_grouped) %>% add_tally(name="n_spots") %>%
               mutate(n_spots=ifelse(problem_areas_grouped=="none", NA, n_spots)), 
         aes(y=sum_umi, x=problem_areas_grouped, fill=n_spots))+
  geom_violin(draw_quantiles= c(.25,.5,.75))+geom_hline(aes(yintercept=100), color="tomato", linewidth=3, alpha=.5)+
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(n=6,"BrBG")[6:1])+
  geom_text(data=filter(cdata, in_tissue==TRUE, problem_areas_grouped!="none") %>% group_by(problem_areas_grouped) %>% tally(name="n_spots"),
            aes(y=10000, label=paste(n_spots,"spots")))+
  geom_text(data=filter(cdata, in_tissue==TRUE, problem_areas_grouped!="none") %>% distinct(problem_areas_grouped, problem_areas_genes.id) %>% 
                group_by(problem_areas_grouped) %>% tally(name="n_areas"),
              aes(y=5000, label=paste(n_areas,"areas"), fill=NULL))+
  geom_text(data=filter(cdata, in_tissue==TRUE, problem_areas_grouped=="none") %>% group_by(problem_areas_grouped) %>% tally(name="n_spots"),
            aes(y=80000, label=paste(n_spots, "spots"), fill=NULL))+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  labs(x="problem area size",y="per-spot UMI (log10 scale)", fill="# spots\nin group",
       title="Library size in problem areas", subtitle="bars in violin indicate Q1, median, and Q3")+
  theme_bw()+theme(axis.title.y=element_text(margin=margin(r=6, unit="pt")), 
                   text=element_text(size=10))

t1 = filter(cdata, in_tissue==TRUE, problem_areas_grouped!="none") %>% group_by(problem_areas_grouped, problem_areas_genes.id) %>% 
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots)
t2 = filter(cdata, in_tissue==TRUE, problem_areas_grouped=="edge") %>% group_by(sample_id, problem_areas_grouped) %>% 
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>% ungroup() 
colnames(t2)[1] = "problem_areas_genes.id"
t3 = bind_rows(filter(t1, problem_areas_grouped!="edge"), t2) %>%
  mutate(n_spots_bin=cut(n_spots, breaks=c(0,20,50,100,500,1100)))
p3 <- ggplot(t3, aes(x=problem_areas_grouped, y=prop_lowumi))+
  geom_hline(aes(yintercept=.5), color="tomato", linewidth=3, alpha=.5)+
  #ggbeeswarm::geom_quasirandom(aes(color=log10(n_spots)))+
  #scale_color_gradient(low="grey", high="#56B1F7")+
  ggbeeswarm::geom_quasirandom(aes(color=n_spots_bin), size=2)+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(n=5, "BrBG")[5:4], "grey", RColorBrewer::brewer.pal(n=5, "BrBG")[2:1]))+
  #scale_color_brewer(palette="BrBG", direction = -1)+
  geom_boxplot(outliers = FALSE, fill=NA)+
  #scale_color_gradientn(colors=c("grey","lightgreen","#56B1F7","#56B1F7"), values=c(0,.1,.5,1))+
  labs(x="problem area size", y="prop. spots in area with <=100 UMI each", color="# spots\nin area",
       title="Problem areas with highly compromised spots", subtitle="each point is a single continuous area")+
  theme_bw()+theme(#plot.margin = margin(t=.5, r=1, b=.5, l=.5, unit="cm"), axis.title.y=element_text(margin=margin(r=10, unit="pt")),
                   text=element_text(size=12))





cdata$problem_areas_grouped2 = ifelse(cdata$problem_areas_genes.size<=20, "small", "flag")
cdata$problem_areas_grouped2 = ifelse(cdata$problem_areas_genes.size==0, "none", cdata$problem_areas_grouped2)
cdata[cdata$problem_areas_binary,"problem_areas_grouped2"] = "remove"
cdata[cdata$true_edges,"problem_areas_grouped2"] = "edge"

cdata$facet_spots = paste(cdata$round, cdata$sample_id)
cdata$facet_spots = ifelse(cdata$brnum=="Br5366", paste(cdata$facet_spots, cdata$brnum), cdata$facet_spots)
x_ordered = sort(unique(cdata$facet_spots))
cdata$x_facets = ""
cdata[cdata$facet_spots %in% x_ordered[1:60],"x_facets"] = "g1"
cdata[cdata$facet_spots %in% x_ordered[61:120],"x_facets"] = "g2"
#colData(spe2)[spe2$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
#colData(spe2)[spe2$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
#colData(spe2)[spe2$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
#colData(spe2)[spe2$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

tmp = filter(cdata, in_tissue==TRUE) %>% group_by(x_facets, facet_spots, problem_areas_grouped2) %>%
  summarise(n_spots=n(), n_lowumi=sum(lowumi))

color.palette2 = c("skyblue","palegoldenrod","lightgreen","black","grey80")
names(color.palette2) <- c("edge","remove","flag","small","none")

ggplot(tmp, aes(x=facet_spots, y=n_lowumi, fill=problem_areas_grouped2))+
  geom_bar(stat="identity", position="stack", color="black", linewidth=.5)+
  scale_y_continuous(limits=c(0,1100), expand=expansion(add=c(0,20)))+
  scale_fill_manual(values=color.palette2)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="# spots", fill="problem\nareas", title="Very low UMI spots across samples")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())

tmp2 = filter(cdata, in_tissue==TRUE) %>% group_by(problem_areas_grouped2, lowumi) %>%
  summarise(n_spots=n()) %>% 
  mutate(problem_areas_grouped2=factor(problem_areas_grouped2, levels=c("edge","remove","flag","small","none")))
ggplot(filter(tmp2, lowumi==TRUE), aes(x="", y=n_spots, fill=problem_areas_grouped2))+
  geom_bar(stat="identity", position="fill")+coord_polar(theta = "y")+
  scale_fill_manual(values=color.palette2)+
  #facet_wrap(vars(lowumi), ncol=2, scales="free_x")+
  theme_void()+labs(title="Very low UMI spots identified for removal", fill="problem\narea group")

ggplot(tmp, aes(x=facet_spots, y=n_lowumi, fill=problem_areas_grouped))+
  geom_bar(stat="identity", position="stack", color="black", linewidth=.5)+
  scale_y_continuous(limits=c(0,1100), expand=expansion(add=c(0,20)))+
  scale_fill_manual(values=cp2)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="# spots", fill="problem\nareas", title="Very low umi problem area spots across samples")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())


