library(dplyr)
library(ggplot2)
set.seed(123)

svg.list = list.files("processed-data/04_feature_selection/per-slide_svgs")
length(svg.list)

svg.df = do.call(rbind, lapply(svg.list, function(x) {
  dummy_slide = strsplit(x, "_")[[1]][1]
  tmp = read.csv(paste0("processed-data/04_feature_selection/per-slide_svgs/",x), row.names=1)
  rownames(tmp) <- NULL
  tmp$dummy_slide= dummy_slide
  tmp
}))

svg.df2 = group_by(svg.df, dummy_slide) %>% mutate(ll=abs(floor(min(LR_stat)))) %>% ungroup() %>%
  mutate(LR_stat_lg10=log10(LR_stat+ll),
         sig=padj<.05)

#################### nnSVG evaluation plots
# n sig genes per slide
p1 <- ggplot(group_by(svg.df2, dummy_slide, sig) %>% tally(), 
             aes(x=dummy_slide, y=n, fill=sig))+
  geom_bar(stat="identity", position="stack")+
  scale_y_continuous("# genes", expand=c(0,0))+
  scale_fill_manual("padj<.05", values=c("grey50","red3"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

# that sig is similarly evaluated
p2 <- ggplot(svg.df2, aes(x=rank, y=LR_stat_lg10, color=sig))+
  geom_point(size=.3)+scale_color_manual("padj<.05", values=c("grey50","red3"))+
  scale_x_continuous(labels=c("0","2k","4k","6k"))+labs(y="log10(LR_stat + c)")+
  facet_wrap(vars(dummy_slide))+theme_bw()

# look at variance 
p3 <- ggplot(svg.df2, aes(x=tau.sq, y=sigma.sq, color=sig))+
  geom_point(size=.3)+scale_color_manual("padj<.05", values=c("grey50","red3"))+
  ggrepel::geom_text_repel(data=filter(svg.df2, sigma.sq>.5 & dummy_slide!="V13B23-333"), 
                           aes(label=gene_name), 
                           size=3, max.overlaps=Inf, min.segment.length = 0,
                           segment.color="grey", color="black")+
  ggrepel::geom_text_repel(data=filter(svg.df2, sigma.sq>.5 & dummy_slide=="V13B23-333"), 
                           aes(label=gene_name), 
                           size=2, max.overlaps=Inf, min.segment.length = 0,
                           segment.color="grey", color="black")+
  labs(x="tau sq. (non-spatial var)", y="sigma sq. (spatial var)")+
  facet_wrap(vars(dummy_slide))+theme_bw()

laymat = matrix(c(1,2,2,3,3), ncol=1)
pdf(file="plots/04_feature_selection/nnSVG_evaluation_plots.pdf", width=7, height=15)
gridExtra::grid.arrange(p1, p2, p3, layout_matrix=laymat)
dev.off()
