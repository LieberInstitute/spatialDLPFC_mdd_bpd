library(ggplot2)
library(dplyr)

cpList = readRDS("plots/colorPalettes.rds")

out1 = read.csv("processed-data/09_SCENIC/cdata_revised-qc-metrics_for-AUCell.csv", row.names=1) %>%
  mutate(smoothed_k9_1663= factor(smoothed_k9_1663, levels=c("L1","L2","L3.4","L5","L6","WM","Vasc","GABA")),
         seurat_label= factor(seurat_label, levels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo")))


ggplot(out1, aes(x=detected, color=seurat_label))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(values=cpList$transfer.bright)+
  coord_cartesian(xlim=c(0,5000))

ggplot(out1, aes(x=detected, color=smoothed_k9_1663))+
  stat_ecdf(linewidth=2)+
  scale_color_manual(values=c(cpList$smoothed.bright, "Vasc"=cpList$transfer.bright[["Micro.Vasc"]], "GABA"=cpList$transfer.bright[["Inhb"]]))+
  coord_cartesian(xlim=c(0,5000))


table(out1$smoothed_k9_1663)
#L1     L2   L3.4     L5     L6     WM   Vasc   GABA 
#64399 100171 143818  72921  88054  39009   4487    341

sapply(levels(out1$seurat_label), function(x) quantile(filter(as.data.frame(out1), seurat_label==x)$detected, probs=c(.01,.05,.1,.5,1)))
#     Micro.Vasc Astro L2.3   L4   Inhb   L5   L6  Oligo
#1%          119   187  406  214  172.0  411  219  158.0
#5%          184   299  757  430  330.1  769  388  289.1
#10%         235   387  978  609  459.0 1005  535  380.2
#50%         521  1067 1980 1419 1312.0 1955 1473  933.0
#100%       5514  8340 8012 7300 7659.0 8585 7088 6151.0


table(out1$seurat_label)
#Micro.Vasc      Astro       L2.3         L4       Inhb         L5         L6 
#     21265     113974      99403      38022      49663      79653      78237 
#Oligo 
#32983 

sapply(levels(out1$smoothed_k9_1663), function(x) quantile(filter(as.data.frame(out1), smoothed_k9_1663==x)$detected, probs=c(.01,.05,.1,.5,1)))
#       L1   L2 L3.4   L5   L6      WM    Vasc   GABA
#1%    128  340  316  430  232  202.08  116.86  526.4
#5%    194  586  583  776  435  340.00  181.30  706.0
#10%   243  748  787 1002  591  434.00  239.60  851.0
#50%   511 1602 1723 1940 1490 1007.00  550.00 1565.0
#100% 4544 8025 8340 8585 6942 5648.00 4272.00 4756.0

