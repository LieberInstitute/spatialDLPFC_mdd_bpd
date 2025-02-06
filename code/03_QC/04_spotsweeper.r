setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(dplyr)
	library(spdep)
	library(scater)
	library(SpotSweeper)
	library(ggplot2)
})
set.seed(123)

#load spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Dim spe:",dim(spe),"\n")

#load coldata and select problem areas
cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
#all edges are appropriately identified except for from the following samples
### V13B23-302_A1 and V13B23-302_C1
cdata$true_edges = ifelse(cdata$slide=="V13B23-302", FALSE, cdata$edge_outlier_genes)

tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>5, !is.na(problem_areas_genes.id))
remove.areas = unique(filter(tmp, prop_lowumi>=.5)$problem_areas_genes.id)
cat("Problem areas to remove based on >=50% spots with UMI<=100:", length(unique(remove.areas)),"\n")


#add edge info and subset
stopifnot(identical(rownames(cdata),rownames(colData(spe))))

spe$edge_outlier_genes = cdata$edge_outlier_genes
spe$true_edges = cdata$true_edges

spe$problem_areas_genes.id = cdata$problem_areas_genes.id
spe$problem_areas_genes.size = cdata$problem_areas_genes.size
spe$problem_areas_binary = spe$problem_areas_genes.id %in% remove.areas

spe$lowumi = spe$sum_umi<=100
spe$remove_problem.areas = spe$problem_areas_binary | spe$true_edges | spe$lowumi

#copy removal vars to cdata for saving after spotsweeper
cdata$problem_areas_binary = spe$problem_areas_binary
cdata$lowumi = spe$lowumi
cdata$remove_problem.areas = spe$remove_problem.areas
spe = spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")

cat("Criteria for spot removal\n")
table(colData(spe)[,c("problem_areas_binary","lowumi","true_edges")])

cat("Total spots for removal:")
table(spe$remove_problem.areas)

spe = spe[,spe$remove_problem.areas==FALSE]
cat("Dim spe (after filtering):",dim(spe),"\n")

#run spotsweeper
cat(format(Sys.time()), "Calculate local outliers (umi counts)...","\n")
spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "umi_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (n genes)...","\n")
spe <- localOutliers(spe, metric = "sum_gene", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "genes_local.outlier"

cat(format(Sys.time()), "Calculate local outliers (mito %)...","\n")
spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "chrM.ratio_local.outlier"

#update coldata
true.spots = cdata$remove_problem.areas==FALSE & cdata$in_tissue==TRUE
stopifnot(identical(rownames(cdata)[true.spots], rownames(colData(spe))))
cdata$umi_local.outlier = NA
cdata[true.spots, "umi_local.outlier"] = colData(spe)$umi_local.outlier
cdata$genes_local.outlier = NA
cdata[true.spots, "genes_local.outlier"] = colData(spe)$genes_local.outlier
cdata$chrM.ratio_local.outlier = NA
cdata[true.spots, "chrM.ratio_local.outlier"] =	colData(spe)$chrM.ratio_local.outlier

#save coldata
write.csv(cdata, "processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=T)
cat("\nUpdated colData saved to: processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv\n")
#save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))


#plot results
#cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=1)

cdata$facet_spots = paste(cdata$round, cdata$sample_id)
cdata$facet_spots = ifelse(cdata$brnum=="Br5366", paste(cdata$facet_spots, cdata$brnum), cdata$facet_spots)

x_ordered = sort(unique(cdata$facet_spots))
cdata$x_facets = ""
cdata[cdata$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
cdata[cdata$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
cdata[cdata$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
cdata[cdata$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

cdata2 = group_by(cdata, x_facets, facet_spots) %>%
        summarise(n_umi=sum(umi_local.outlier, na.rm=T), n_genes=sum(genes_local.outlier, na.rm=T),
                n_chrM.ratio=sum(chrM.ratio_local.outlier, na.rm=T)) %>%
        tidyr::pivot_longer(c("n_umi","n_genes","n_chrM.ratio"),
                names_to="outlier_type", names_prefix="n_", values_to="n_spots")

#also compute the total number of spots excluded (since some spots were flagged by multiple thresholds
cdata3 = group_by(cdata, x_facets, facet_spots) %>%
        summarise(any.outlier=sum(umi_local.outlier | genes_local.outlier | chrM.ratio_local.outlier, na.rm=T)) %>%
        tidyr::pivot_longer("any.outlier", names_to="outlier_type", values_to="n_spots")

cdata4 = bind_rows(cdata2, cdata3) %>% mutate(is_total = outlier_type=="any.outlier")

p1 <- ggplot(cdata4, aes(x=facet_spots, y=n_spots))+
        geom_bar(data=filter(cdata4, is_total==FALSE), aes(fill=outlier_type),
                stat="identity", position="stack", color="black", alpha=.5, linewidth=.5)+
        geom_bar(data=filter(cdata4, is_total==TRUE), stat="identity", fill="grey30", color="grey30", width=.5)+
        geom_text(data=filter(cdata4, is_total==TRUE), aes(y=0, label=n_spots), vjust=0, color="white", size=3)+
        #scale_y_continuous(expand=expansion(add=c(0,20)))+
        #scale_fill_manual(values=color.palette2)+
        facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
        labs(x="", y="# spots", fill="QC", title="Spotsweeper outliers")+theme_bw()+
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
                strip.placement = "inside", strip.text=element_blank(),
                strip.background = element_blank())

ggsave(filename="plots/03_QC/spotsweeper_outliers_barplot.png", p1, bg="white", units="in", height=12, width=9)
cat("\nPlot saved to: plots/03_QC/spotsweeper_outliers_barplot.png\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
