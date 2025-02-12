library(raster)
library(dplyr)
library(ggplot2)
source("code/03_QC/edgeDetection_finalized/raster_edge_functions.r")

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)

#individual sample walk thru for plotting
#test = cdata[cdata$sample_id=="V13B23-342_D1",c("array_row","array_col","genes_3MAD.outlier_binary")]
#test2 = cdata[cdata$sample_id=="V13B23-342_D1",c("array_row","array_col","in_tissue")]
test = cdata[cdata$sample_id=="V13Y10-023_A1",c("array_row","array_col","genes_3MAD.outlier_binary")]
test2 = cdata[cdata$sample_id=="V13Y10-023_A1",c("array_row","array_col","in_tissue")]
colnames(test2) <- c("x","y","in_tissue")
odds = seq(1,max(test[,"array_col"]), by=2)

t1 = rasterFromXYZ(test)

cp1 = c("grey","black","white")
names(cp1) = c("0","1","off")
t1.df = as.data.frame(t1, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
t1.df = merge(t1.df, test2)
ggplot(mutate(t1.df, fill1=ifelse(in_tissue==FALSE, "off", as.character(genes_3MAD.outlier_binary))), 
       aes(x=x, y=y, fill=fill1))+
  geom_tile()+scale_fill_manual(values=cp1)+theme_void()

p1o <- ggplot(filter(t1.df, !is.na(genes_3MAD.outlier_binary)) %>% 
         mutate(fill1=ifelse(in_tissue==FALSE,"off",as.character(genes_3MAD.outlier_binary))),
       aes(x=x, y=y, fill=fill1))+
  geom_tile()+scale_fill_manual("", values=cp1)+
  theme_void()+ggtitle("raw orig")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)
p1s <- ggplot(filter(t1.df, !is.na(genes_3MAD.outlier_binary)) %>% 
                mutate(fill1=ifelse(in_tissue==FALSE,"off",as.character(genes_3MAD.outlier_binary))),
              aes(x=x, y=y2, fill=fill1))+
  geom_tile(color="white", linewidth=.3)+scale_fill_manual("", values=cp1)+
  theme_void()+ggtitle("raw shifted")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)


t2 = focal_transformations(t1)
t2.df = as.data.frame(t2, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
t2.df = merge(t2.df, test2)
p2o <- ggplot(mutate(t2.df, fill1=ifelse(in_tissue==FALSE, "off", as.character(layer))), 
       aes(x=x, y=y, fill=fill1))+
  geom_tile()+scale_fill_manual(values=cp1)+
  theme_void()+ggtitle("transformed orig")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)
p2s <- ggplot(mutate(t2.df, fill1=ifelse(in_tissue==FALSE, "off", as.character(layer))), 
              aes(x=x, y=y2, fill=fill1))+
  geom_tile(color="white", linewidth=.3)+scale_fill_manual(values=cp1)+
  theme_void()+ggtitle("transformed shifted")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)



c1 = clump(t2, direction=8)
c1.df = as.data.frame(c1, xy=T) %>% mutate(y2 = ifelse(y %in% odds, y-1, y))
c1.df = merge(c1.df, test2)
#recode clump is so that they are more randomized
set.seed(123)
rnum = sample(1:max(c1.df$clumps, na.rm=T), max(c1.df$clumps, na.rm=T))
c1.df$clumps2 = factor(c1.df$clumps, levels=1:max(c1.df$clumps, na.rm=T), 
                       labels=paste(rnum, 1:max(c1.df$clumps, na.rm=T), sep="_"))
p3o <- ggplot(mutate(c1.df, fill1=ifelse(is.na(clumps2), "0", as.character(clumps2))) %>%
           mutate(fill1=ifelse(in_tissue==FALSE,"off",fill1)), 
       aes(x=x, y=y, fill=fill1))+
  geom_tile()+scale_fill_manual(values=c("grey",viridisLite::turbo(n=max(c1.df$clumps, na.rm=T)),"white"))+
    #scale_fill_manual(values=c(viridisLite::turbo(n=max(c1.df$clumps, na.rm=T)),"white"), na.value="grey")+
  theme_void()+ggtitle("clumped orig")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)

p3s <- ggplot(mutate(c1.df, fill1=ifelse(is.na(clumps2), "0", as.character(clumps2))) %>%
         mutate(fill1=ifelse(in_tissue==FALSE,"off",fill1)), 
       aes(x=x, y=y2, fill=fill1))+
  geom_tile(color="white", linewidth=.3)+scale_fill_manual(values=c("grey",viridisLite::turbo(n=max(c1.df$clumps, na.rm=T)),"white"))+
  theme_void()+ggtitle("clumped shifted")+
  theme(legend.position="none", plot.title=element_text(margin=margin(t=.5, l=.5,unit="cm")),
        aspect.ratio=1)


gridExtra::grid.arrange(p1o, p2o, p3o, p1s, p2s, p3s, ncol=3)
