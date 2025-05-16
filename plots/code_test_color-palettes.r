## COLOR PALETTE DESIGN HELP
#
## How to get equal luminance: https://medium.com/sketch-app-sources/mixing-colours-of-equal-luminance-part-2-3e10c07c947c
#
## HSL adjustment of multi color palette: https://supercolorpalette.com/
#  
## I like these sliders for adjusting 1 color: https://www.color-hex.com/color/b1e2b5
#### gives you in HSV so I use this to convert to HEX: https://atmos.style/color-converter/hex-to-lch
#### this one also converts to CYMK which is what is needed for printing: https://convertacolor.com/ 
#
## I like this site for exploring HEX variants of a single color: https://www.color-hex.com/color/f5d29e
#
## Color blindness check: https://www.color-blindness.com/coblis-color-blindness-simulator/
#

dx.pal = c("grey50","#9e771b","#1b9e77")#oops very not colorblind friendly
#try taking some values from viridis gradients
#viridis::magma(n=5)
#"#000004FF" "#51127CFF" "#B63679FF" "#FB8861FF" "#FCFDBFFF"
dx.pal = c("grey50","#FB8861FF","#51127CFF")
#modify values for brightness
dx.pal = c("grey50","#FB8861FF","#7e3fa9")
dx.pal = c("#7f7f7f","#FB8861FF","#9260b2") #sRGB friendly
names(dx.pal) = c("NTC","MDD","BPD")


### easter color palette
easter.pal1 = c("#FF8A97","#F5D29E","#B1E2B5","#B1BEEC","#D7A5EE","#F5EE9E","#E2B1C6","grey90")
names(easter.pal1) = c("Vasc","L1","L2","L3","GABA","L5","L6","WM")
easter.pal2 = c("#ff3e53","#eeb259","#79cd80","#738bdd","#B58AFF","#eee259","#cd799d","grey")
names(easter.pal2) = c("Vasc","L1","L2","L3","GABA","L5","L6","WM")



###earthy color palette
earthy.pal1 = c("#d05e46","#F5D29E","#b3c7ac","#a3c7e4","#b9a7c9","#e5e8a0","#ef9e9f","#ede8e0")
names(earthy.pal1) = c("Vasc","L1","L2","L3","GABA","L5","L6","Whsl(207,55%,77%)M")
earthy.pal2 = c("#911223","#cfa45c","#5D9940","#5095CD","#9377AC","#ddc94e","#E45C5F","#D1C4B0")
# "#699e59" -- swap out green for
names(earthy.pal2) = c("Vasc","L1","L2","L3","GABA","L5","L6","WM")

saveRDS(list(dx.pal= dx.pal, earthy.pal1= earthy.pal1, earthy.pal2= earthy.pal2),
        "plots/colorPalettes.rds")

#fifth element
element.list = list(c("#2C87A1","#5D9940","#E34611","#A76E51","#DBDA2F","#4EBBC7","#85A0A0"),
                    c("#3495CA","#68BA53","#F27D35","#BD936F","#D1E457","#71BDD5","#9EB3B6"),
                    c("#599AD7","#7ACB75","#F7AF62","#CEB590","#D0EC81","#96C6E1","#BAC8CD"),
                    c("#80A5E2","#98DA9C","#FAD590","#DFD2B1","#D7F4AC","#BCD5ED","#D6DEE3"))


#palettes based on images in my library
#froggey
frog.pal <- c("#2257A7","#8E7360","#43422D","#8657BF","#E4775D","#F2CF74","#7BAAEE","#D2C9BA")


#geology pride flag
geo.pal <- c("#7F8D11","#941D61","#C13E20","#475FA5","#A0C255","#EEBC4A","#E872AA","#93D3F6")


#palm springs
palm.pal <- c("#4E1304","#004871","#4C3B47","#CD624C","#716DAA","#AB9C75","#FFD351","#93D3F6")


#starry night
star.pal <- c("#222533","#546129","#9E4D12","#496E9A","#A38A6E","#D3D051","#9EB5BF","#ECE7DD")


#spring flowers
spring.pal <- c("#463929","#A04571","#6C9066","#90D542","#E077B6","#B0C1B4","#F8F7A9","#E2E5F5")




library(SpatialExperiment)
library(HDF5Array)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(escheR)

set.seed(123)

#color.palettes = list("light"= easter.pal1, "bright"= easter.pal2)
color.palettes = list("light"= earthy.pal1, "bright"= earthy.pal2)
#test color palettes with:

# 4 spot plot of just domains
# spot plot of gradient (high expr)
# spot plot of gradient (med expr)
# spot plot of gradient (low expr)
# stacked/filled bar plot domain (light and bright)
# boxplot (domain) / violin (domain)
# box plot (dx)
# heatmap greyscale
# heatmap diverging
# domain color gradient?


#spot plots
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
spotdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
spotdata = spotdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(spotdata)))

spotdata = spotdata[spotdata$precast_k9_1663_f!="low UMI",]
spe = spe[,rownames(spotdata)]
spe$precast_k9_1663 = factor(spotdata$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
table(spe$precast_k9_1663, useNA="ifany")

#domain spot plots
sub.samples = c("329-A1", "329-C1", "309-C1", "342-D1")
sub.samples = paste0("V13B23-", gsub("-","_", sub.samples))

spe_sub = spe[,spe$sample_id %in% sub.samples]

mod_spatialCoords = spatialCoords(spe_sub)
for (i in unique(spe_sub$sample_id)) {
  tmp = mod_spatialCoords[colData(spe_sub)$sample_id==i,]
  mod_spatialCoords[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
  mod_spatialCoords[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
}

p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords[,1], y_coord=mod_spatialCoords[,2],
                sample_id="sample_id", annotate="precast_k9_1663", point_size=.1)+
  scale_color_manual("PRECAST\ncluster", values=color.palettes[[2]])+
  facet_wrap(vars(sample_id), ncol=4)+
  guides(colour = guide_legend(nrow= 1, override.aes = list(size = 3)))+
  theme(strip.text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), 
        panel.border=element_rect(fill=NA, color=NA),
        legend.position="bottom", legend.key.size = unit(6,"pt"),
        #legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"))#,
        #legend.box.margin = margin(0,2,0,2,"pt"))

#gex spot plots
spe_sub2 = spe[,spe$sample_id=="V13B23-329_A1"]

spot.genes = c("PCDH8","RORB","PCP4","CLSTN2",#general gradients
               "DGKA","TNNT2","TRABD2A","OPRK1")#specific 
for(i in spot.genes) {
  spe_sub2[[i]] = logcounts(spe_sub2)[rowData(spe_sub2)$gene_name==i,]
}

plist <- lapply(c("PCDH8","RORB","TNNT2","TRABD2A"), function(x) {
  p = make_escheR(spe_sub2) %>%
    add_ground(var="precast_k9_1663", stroke=.2, point_size=.5) %>% 
    add_fill(var=x, point_size = .5)
  p+scale_color_manual("", values=color.palettes[[1]])+
    scale_fill_gradient("log2\nCPM", low="white",high="black")+
    labs(title=x)+theme(text=element_text(size=10), plot.title=element_text(face="italic"),
                        legend.position="none")
})



#other plots
load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")
cdata = as.data.frame(colData(spe_pseudo))
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))


p2 <- ggplot(group_by(cdata, cond_sex, precast_k9_1663) %>% summarise(n_total=sum(nspots)), 
             aes(x=cond_sex, y=n_total, fill=precast_k9_1663))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_fill_manual(values=color.palettes[[2]])+
  guides(fill = guide_legend(nrow= 1))+
  labs(y="# spots", x="", fill="PRECAST\ncluster", title="Cluster representation")+
  theme_bw()+theme(legend.position="none", axis.title.x=element_blank())

p3 <- plotExpression(spe_pseudo, features=c("SLC32A1"), colour_by = "precast_k9_1663",
               x="precast_k9_1663", swap_rownames = "gene_name", assay.type = "logcounts")+
  scale_color_manual(values=color.palettes[[2]])+labs(y="log2 CPM")+
  theme(legend.position="none", axis.title.x=element_blank())


#rasd1 condition plot
cdata$rasd1 = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name=="RASD1",]
p4 <- ggplot(cdata, aes(x=condition, y=rasd1, fill=condition))+
  geom_boxplot(color="black", outlier.size=.5, linewidth=.5)+
  scale_fill_manual(values=dx.pal, guide="none")+
  labs(y="log2 CPM", title="RASD1")+
  theme_bw()+theme(text=element_text(size=9), plot.title=element_text(face="italic"),
                   panel.grid.minor = element_blank())

#heatmap
marker.genes = c("CLDN5", "FABP7", "HPCAL1", "KCNH5", "RORB", "SLC32A1", "PCP4", "CLSTN2", "MOXD1", "SHTN1")

hmp = plotGroupedHeatmap(spe_pseudo, features=marker.genes, group="precast_k9_1663", 
                         swap_rownames="gene_name", 
                         scale = T, center=T, cluster_cols=F, cluster_rows=F,
                         angle_col=0, border_color="white",
                         colour = RColorBrewer::brewer.pal(7, "Greys"), main="Cluster markers")
                         #colour = RColorBrewer::brewer.pal(7, "RdGy")[7:1])


lay.mat = rbind(c(1,1,1,2),
                c(3,3,4,5),
                c(3,3,6,7),
                c(8,8,9,9))

ggsave("plots/earthy_color_palette.png", #"plots/earthy_color_palette.png", 
       gridExtra::grid.arrange(p1, p4,
                               hmp[[4]], 
                               plist[[1]], plist[[2]], plist[[3]], plist[[4]], 
                               #plist[[5]], plist[[6]], plist[[7]], plist[[8]],
                               p2, p3,
                               layout_matrix=lay.mat),
       bg="white", width=8, height=9)
