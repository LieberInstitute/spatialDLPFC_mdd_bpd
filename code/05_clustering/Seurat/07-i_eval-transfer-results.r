library(Seurat)
set.seed(123)

##################
### found that the kweights=50 was cleaner (less random spotting) than kweights=20
##################

### PCs = 20
##### kweight=20
res1 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20_red-precast-kweight-20-low-res.csv", row.names=1)
hist(res1$prediction.score.max, breaks=50) #pretty damn good
table(res1$predicted.id)
#Astro       Inhb         L2         L3         L4         L5         L6 
#79827      37106      29659      51848     109160      87080      65290 
#Micro/Vasc      Oligo 
#.....30714      44564 
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
identical(rownames(res1), rownames(cdata))

cdata$predicted.id = res1$predicted.id
(t1 = table(cdata[,c("precast_k9_1663_f","predicted.id")]))
#..................predicted.id
#precast_k9_1663_f Astro  Inhb    L2    L3    L4    L5    L6 Micro/Vasc Oligo
#..........GABA        0  7002     0     0     1     2     0          0     0
#..........L1      28183 12494   158     2  1400    54  2089      13894   352
#..........L2      21451  9286 27631 23998 10032   381  2704       1080     0
#..........L3      10638  3148   477 27790 73224 18334  6205        164    11
#..........L5       2783   377   469     0  6337 60540   362         80     1
#..........L6       7741  1183    42    35 12710  6035 50780       2997  2355
#..........low UMI  1011  3365   857    18  5323  1043  1543       3328  9761
#..........Vasc      582   231     8     1    28   569  1501       8963    69
#..........WM       7438    20    17     4   105   122   106        208 32015
round(t1/rowSums(t1),2)*100
#..................predicted.id
#precast_k9_1663_f Astro Inhb  L2  L3  L4  L5  L6 Micro/Vasc Oligo
#..........GABA        0  100   0   0   0   0   0          0     0
#..........L1         48   21   0   0   2   0   4         24     1
#..........L2         22   10  29  25  10   0   3          1     0
#..........L3          8    2   0  20  52  13   4          0     0
#..........L5          4    1   1   0   9  85   1          0     0
#..........L6          9    1   0   0  15   7  61          4     3
#..........low UMI     4   13   3   0  20   4   6         13    37
#..........Vasc        5    2   0   0   0   5  13         75     1
#..........WM         19    0   0   0   0   0   0          1    80


##### kweight=50
res2 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20_red-precast-kweight-50-low-res.csv", row.names=1)
hist(res2$prediction.score.max, breaks=50)
table(res2$predicted.id)
#Astro       Inhb         L2         L3         L4         L5         L6 
#96305      36948      26306      52024     103361      84841      66665 
#Micro/Vasc      Oligo 
#.....23704      45094
cdata$predicted.id2 = res2$predicted.id
(t2 = table(cdata[,c("precast_k9_1663_f","predicted.id2")], useNA="ifany"))
#predicted.id2
#precast_k9_1663_f Astro  Inhb    L2    L3    L4    L5    L6 Micro/Vasc Oligo
#..........GABA        0  7004     0     0     0     1     0          0     0
#..........L1      33178 12615     7     0   513    74   965      11021   253
#..........L2      26611  9883 25699 24032  8507    23  1681        127     0
#..........L3      14901  2480    13 27983 77857 14378  2379          0     0
#..........L5       2732   232    18     0  2767 65144    55          1     0
#..........L6      10596   961     3     7  7664  4329 59087        523   708
#..........low UMI   929  3396   564     0  5959   656  1043       2583 11119
#..........Vasc      722   357     1     2    33   177  1423       9171    66
#..........WM       6636    20     1     0    61    59    32        278 32948
round(t2/rowSums(t2),2)*100
#predicted.id2
#precast_k9_1663_f Astro Inhb  L2  L3  L4  L5  L6 Micro/Vasc Oligo
#..........GABA        0  100   0   0   0   0   0          0     0
#..........L1         57   22   0   0   1   0   2         19     0
#..........L2         28   10  27  25   9   0   2          0     0
#..........L3         11    2   0  20  56  10   2          0     0
#..........L5          4    0   0   0   4  92   0          0     0
#..........L6         13    1   0   0   9   5  70          1     1
#..........low UMI     4   13   2   0  23   2   4         10    42
#..........Vasc        6    3   0   0   0   1  12         77     1
#..........WM         17    0   0   0   0   0   0          1    82





### PCs = 30
##### kweight=20
res3 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-20-low-res.csv", row.names=1)
hist(res3$prediction.score.max, breaks=50) #pretty damn good
table(res3$predicted.id)
#Astro       Inhb         L2         L3         L4         L5         L6 
#97405      52116      14192      82368      50057      81756      85900 
#Micro/Vasc      Oligo 
#.....29902      41552 
identical(rownames(cdata), rownames(res3))
cdata$predicted.id3 = res3$predicted.id
  #ifelse(res3$prediction.score.max<.5, NA, res1$predicted.id)
(t1 = table(cdata[,c("precast_k9_1663_f","predicted.id3")], useNA="ifany"))
#predicted.id
#precast_k9_1663_f Astro  Inhb    L2    L3    L4    L5    L6 Micro/Vasc Oligo
#..........GABA       10  6976     0     0     0    19     0          0     0
#..........L1      30024  9841   117     1   185    39  7393      10902   124
#..........L2      29440 13621 11744 33773  3814    30  2736       1405     0
#..........L3      15239 12559  1502 45526 36800 16674 10938        293   460
#..........L5       3545  2068    23     1  3265 59359  1956        348   384
#..........L6       8713  3261    41  2881  2477  3611 58821       2513  1560
#..........low UMI  1468  3509   765   106  3513  1105  3117       4141  8525
#..........Vasc      842   230     0    76     0   785   477       9388   154
#..........WM       8124    51     0     4     3   134   462        912 30345
round(t1/rowSums(t1),2)*100
#predicted.id
#precast_k9_1663_f Astro Inhb  L2  L3  L4  L5  L6 Micro/Vasc Oligo
#..........GABA        0  100   0   0   0   0   0          0     0
#..........L1         51   17   0   0   0   0  13         19     0
#..........L2         30   14  12  35   4   0   3          1     0
#..........L3         11    9   1  33  26  12   8          0     0
#..........L5          5    3   0   0   5  84   3          0     1
#..........L6         10    4   0   3   3   4  70          3     2
#..........low UMI     6   13   3   0  13   4  12         16    32
#..........Vasc        7    2   0   1   0   7   4         79     1
#..........WM         20    0   0   0   0   0   1          2    76


##### kweight=50
res4 <- read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
hist(res4$prediction.score.max, breaks=50)
table(res4$predicted.id)
#.Astro       Inhb         L2         L3         L4         L5         L6 
#115199      52740       8770      90981      40983      80303      80688 
#Micro/Vasc      Oligo 
#.....24519      41065 
#### ALMOST NO L2

identical(rownames(cdata), rownames(res4))
cdata$predicted.id4 = res4$predicted.id

(t4 = table(cdata[,c("precast_k9_1663_f","predicted.id4")], useNA="ifany"))
#..................predicted.id4
#precast_k9_1663_f Astro  Inhb    L2    L3    L4    L5    L6 Micro/Vasc Oligo
#..........GABA        0  7004     0     0     0     1     0          0     0
#..........L1      34058  8252     9     0    34     9  7060       9187    17
#..........L2      34094 15121  8113 35050  2834     4  1254         93     0
#..........L3      22157 13173   257 54777 31867 12419  5161        131    49
#..........L5       3031  1399     0     2   897 64686   810         25    99
#..........L6      11405  3610     2  1083  1590  1866 63133        801   388
#..........low UMI  1433  3924   389    21  3759   838  2737       3798  9350
#..........Vasc      920   208     0    48     0   435   315       9958    68
#..........WM       8101    49     0     0     2    45   218        526 31094
round(t4/rowSums(t4),2)*100
#predicted.id4
#precast_k9_1663_f Astro Inhb  L2  L3  L4  L5  L6 Micro/Vasc Oligo
#..........GABA        0  100   0   0   0   0   0          0     0
#..........L1         58   14   0   0   0   0  12         16     0
#..........L2         35   16   8  36   3   0   1          0     0
#..........L3         16    9   0  39  23   9   4          0     0
#..........L5          4    2   0   0   1  91   1          0     0
#..........L6         14    4   0   1   2   2  75          1     0
#..........low UMI     5   15   1   0  14   3  10         14    36
#..........Vasc        8    2   0   0   0   4   3         83     1
#..........WM         20    0   0   0   0   0   1          1    78


#test 2 pseudobulk datasets
### one from pc20 kweight=50
### one from pc30 kweight=50
### for both pseudobulk datasets, combine L2 and L3 into one
pc20 = factor(res2$predicted.id, levels=c("Micro/Vasc","Astro",
                                          "L2","L3",
                                          "L4","L5","L6",
                                          "Oligo","Inhb"),
              labels=c("Micro/Vasc","Astro","L2/3","L2/3",
                       "L4","L5","L6",
                       "Oligo","Inhb"))
pc30 = factor(res4$predicted.id, levels=c("Micro/Vasc","Astro",
                                          "L2","L3",
                                          "L4","L5","L6",
                                          "Oligo","Inhb"),
              labels=c("Micro/Vasc","Astro","L2/3","L2/3",
                       "L4","L5","L6",
                       "Oligo","Inhb"))
### differences will be in L4, and inhibitory
rbind(table(pc20, useNA="ifany"),
      table(pc30, useNA="ifany"))
#Micro/Vasc  Astro  L2/3     L4    L5    L6 Oligo  Inhb
#.....23704  96305 78330 103361 84841 66665 45094 36948
#.....24519 115199 99751  40983 80303 80688 41065 52740




### spot plots
library(HDF5Array)
library(SpatialExperiment)
library(ggspavis)

cpList = readRDS("plots/colorPalettes.rds")

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
identical(rownames(colData(spe)), rownames(cdata))

spe$precast_k9_1663 = cdata$precast_k9_1663_f
class(spe$precast_k9_1663)
spe$predicted.id = cdata$predicted.id
spe$predicted.id0 = cdata$predicted.id0
spe$predicted.id4 = cdata$predicted.id4
vistoseg.samples = c("308-D1", "329-A1", "342-B1", "332-B1", #NTC M
                     "342-A1", "332-A1", "327-C1", "329-B1",#NTC F
                     "382-C1", "309-D1","352-A1", "329-C1",  #MDD M
                     "382-D1", "352-B1", "309-C1", "380-B1", #"329-D1", #MDD F
                     "382-A1", "352-C1", "342-D1", "327-B1", #BPD M
                     "382-B1", "308-A1", "309-A1", "327-A1") #BPD F

vistoseg.samples = paste0("V13B23-", gsub("-","_", vistoseg.samples))

spe_sub = spe[,spe$sample_id %in% vistoseg.samples]

spe_sub$cond_sex = factor(paste(spe_sub$condition, spe_sub$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))
spe_sub$facet_col = factor(spe_sub$sample_id, levels=vistoseg.samples, labels=paste0("C",rep(c(1:4), 6)))

mod_spatialCoords = spatialCoords(spe_sub)
for (i in unique(spe_sub$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

low.res.pal = c("Astro"="#cfa45c","Micro/Vasc"="#911223",
                "Inhb"="#9377AC",
                "L2"="#5D9940", "L3"="#5095CD", 
                "L4"="#85A0A0",
                "L5"="#ddc94e","L6"="#E45C5F","L6b"="#ef9e9f",
                "Oligo"="#D1C4B0")

low.res.pal2 = c("Astro"="#cfa45c","Micro/Vasc"="#911223",
                "Inhb"="#9377AC",
                "L2"="#375b26", "L3"="#5095CD", 
                "L4"="#c2cfcf",
                "L5"="#ddc94e","L6"="#E45C5F","L6b"="#ef9e9f",
                "Oligo"="grey90")

p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
                sample_id="sample_id", annotate="predicted.id4", point_size=.1)+
  scale_color_manual("predicted\ncluster", values=low.res.pal2)+
  facet_grid(rows=vars(cond_sex), cols=vars(facet_col), switch="y")+
  geom_text(aes(x=0, y= -60, 
                label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13,
                              .data[["sample_id"]], "")), hjust=0, vjust=0, 
            color="black", size=2)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

ggsave("plots/05_clustering/Seurat/VistoSeg-examples_qual-genes-k-50-anchors-pc30_transfer-precast-red-kweight-50-low-res_spot-plots-recolor.png", 
       p1,
       bg="white", width=8, height=10)


#lower quality examples
lowq.samples = c("V13B23-282","V13B23-301","V13B23-302","V13B23-310","V13B23-311","V13B23-403")
spe_sub2 = spe[,spe$slide %in% lowq.samples]

spe_sub2$facet_row = factor(spe_sub2$slide, levels= lowq.samples, labels= gsub("-","\n",lowq.samples))
spe_sub2$facet_col = factor(spe_sub2$array, levels= c("A1","B1","C1","D1"))

mod_spatialCoords = spatialCoords(spe_sub2)
for (i in unique(spe_sub2$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub2)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub2)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub2)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

p2 <- plotSpots(spe_sub2, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
                sample_id="sample_id", annotate="predicted.id4", point_size=.1)+
  scale_color_manual("predicted\ncluster", values=low.res.pal2)+
  facet_grid(rows=vars(facet_row), cols=vars(facet_col), switch="y")+
  #geom_text(aes(x=0, y= -60, 
  #              label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13,
  #                            .data[["sample_id"]], "")), hjust=0, vjust=0, 
  #          color="black", size=2)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        strip.text.x = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

#ggsave("plots/05_clustering/Seurat/low-qual-examples_k-100-anchors-transfer-low-res_spot-plots.png", 
ggsave("plots/05_clustering/Seurat/low-qual-examples_qual-genes-k-50-anchors-pc30_transfer-precast-red-kweight-50-low-res_spot-plots-recolor.png", 
       p2,
       bg="white", width=8, height=10)

