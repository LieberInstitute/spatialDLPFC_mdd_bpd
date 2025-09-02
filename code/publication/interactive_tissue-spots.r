library(HDF5Array)
library(SpatialExperiment)
library(DelayedArray)
library(escheR)
library(ggspavis)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
spe_sub = spe[,spe$sample_id=="V13B23-329_A1"]

p <- plotVisium(spe_sub, spots=FALSE, image=TRUE)

p <- p %>% add_ground(var="in_tissue", stroke=.2, point_size=.8)
p+scale_color_manual(values=c("grey50"), guide="none")

ggsave(file="plots/publication/tissue_spots.png", 
       p+scale_color_manual(values=c("black"), guide="none"),
       bg="white", height=5, width=5)


head(demo)
table(demo[,c("array","cond_sex")])

df = tidyr::pivot_wider(demo[,c("array","cond_sex","slide")], names_from="array", values_from="cond_sex")
df = df[,c("slide","A1","B1","C1","D1")]
group_by(df, A1, B1, C1, D1) %>% tally()
