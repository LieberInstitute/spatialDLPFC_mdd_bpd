library(SpatialExperiment)
library(HDF5Array)
library(ggspavis)
library(ggplot2)
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
cdata2 = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")

identical(rownames(cdata), colnames(spe))

spe$genes_3MAD.outlier_binary = cdata$genes_3MAD.outlier_binary
spe$lg10.umi = log10(spe$sum_umi)
spe$umi100 = spe$sum_umi<=100
spe$true_edges = cdata2$true_edges
spe$problem_areas_binary = cdata2$problem_areas_binary
spe$plot_problem.areas = spe$true_edges | spe$problem_areas_binary

plist = list()

example_sections = c("V13B23-334_C1", "V13B23-328_A1")

for(i in example_sections) {
  
  spe_sub = spe[,spe$sample_id==i]
  spe_sub = spe_sub[,spe_sub$in_tissue==T]
  
  plist[[length(plist)+1]] = plotSpots(spe_sub, annotate="lg10.umi", point_size = .9)+ggtitle("Sum UMI (log10)")+
    theme(panel.background = element_rect(fill="grey50"),
          legend.box.spacing = unit(6,"pt"))
  plist[[length(plist)+1]] = plotSpots(spe_sub, annotate="sum_gene", point_size = .9)+
    scale_color_gradient("", low="#F2F2F2", high="#0000FF", labels=function(x) paste0(x/1000,"k"))+ggtitle("Detected genes")+
    theme(panel.background = element_rect(fill="grey50"), 
          legend.box.spacing = unit(3,"pt"))
  plist[[length(plist)+1]] = plotSpots(spe_sub, annotate="expr_chrM_ratio", point_size = .9)+ggtitle("Mito fraction")+
    theme(panel.background = element_rect(fill="grey50"),
          legend.box.spacing = unit(2,"pt"))
  
  
  spe_sub2 = spe_sub[,spe_sub$plot_problem.areas==F]
  plist[[length(plist)+1]] = plotSpots(spe_sub2, annotate="lg10.umi", point_size = .9)+ggtitle("Sum UMI (log10)")+
    theme(panel.background = element_rect(fill="grey50"),
          legend.box.spacing = unit(2,"pt"))
  plist[[length(plist)+1]] = plotSpots(spe_sub2, annotate="sum_gene", point_size = .9)+
    scale_color_gradient("", low="#F2F2F2", high="#0000FF", labels=function(x) paste0(x/1000,"k"))+ggtitle("Detected genes")+
    theme(panel.background = element_rect(fill="grey50"), 
          legend.box.spacing = unit(3,"pt"))
  plist[[length(plist)+1]] = plotSpots(spe_sub2, annotate="expr_chrM_ratio", point_size = .9)+ggtitle("Mito fraction")+
    theme(panel.background = element_rect(fill="grey50"),
          legend.box.spacing = unit(2,"pt"))
  
  
  
}


ggsave(file="plots/publication/global_outliers.pdf", 
       gridExtra::marrangeGrob(plist, ncol=3, nrow=2, layout_matrix=matrix(seq_len(2 * 3), nrow = 2, ncol = 3, byrow = T),
                               top=quote(example_sections[[g]])),
       bg="white", height=10, width=12)

#qc outcomes
cdata2 = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
cdata3 = cdata2[cdata2$in_tissue,]
cdata3$cond_sex = factor(paste(cdata3$condition, cdata3$sex), levels=c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))
tmp = group_by(cdata3, condition, sex, cond_sex, sample_id, remove_spots) %>% tally() %>% 
  group_by(sample_id) %>% mutate(total_spots=sum(n)) %>%
  ungroup() %>% mutate(prop_spots=n/total_spots, 
                       remove_spots= factor(as.character(remove_spots), levels=c("FALSE","TRUE"),
                                            labels=c("Keep","Remove")))

p3 <- ggplot(tmp, aes(x=remove_spots, y=prop_spots, color=condition, shape=sex))+
  ggbeeswarm::geom_quasirandom()+scale_color_manual(values=cpList$dx.pal)+
  facet_wrap(vars(cond_sex), ncol=6)+
  labs(x="",y="prop. of spots")+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"))

ggsave(file="plots/publication/qc_excluded.png", 
       p3,
       bg="white", height=3, width=8)
