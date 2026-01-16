setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(grid)
	library(ggrastr)
	library(UpSetR)
})
set.seed(123)

#results_set = "smoothed-k9-1663"
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")


comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

covars = c("none","pc3 only","detected","nspots","age","BMI","Smoking","RIN","PMI")
names(covars) <- covars

covars.colors = c("grey","black","#E872AA","#2257A7","#E4775D","#A0C255","#EEBC4A","#A38A6E","#93D3F6")
names(covars.colors) = names(covars)

# layer adjusted
## F test
f1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                     results_set, "_rev-gene-input_F-test.csv"), row.names=1) %>%
  mutate(covar="none")
f2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                     results_set, "_rev-gene-input_F-test.csv")) %>%
  mutate(covar="pc3 only")
f3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                       results_set, "_rev-gene-input_F-test_all-covar-models.csv"))

la.all = bind_rows(f1, f2, f3) %>% 
  mutate(covar=factor(covar, levels=covars))

#filter(la.all, adj.P.Val<.05) %>% group_by(covar) %>% tally()
#filter(la.all, adj.P.Val<.05) %>% group_by(covar) %>% slice_min(F) %>%
#  select(covar, F, P.Value, adj.P.Val)

p1 <- ggplot(filter(la.all, adj.P.Val<.05), aes(x=covar, fill=covar))+
  geom_bar(stat="count", color="black", linewidth=.5)+
  geom_hline(aes(yintercept=nrow(filter(la.all, adj.P.Val<.05, covar=="pc3 only"))), lty=2)+
  scale_fill_manual(values=covars.colors, guide="none")+
  labs(y="# genes", title="Layer-adjusted F test adj p<.05", subtitle=results_set,
       x="covariate in model")+
  theme_minimal()

covars1 = c("nspots","detected","RIN")
names(covars1) <- covars1
checkList1 <- lapply(c("pc3 only", covars1), function(x) filter(la.all, adj.P.Val<.05, covar==x)$gene_id)
names(checkList1) <- c("pc3", covars1)

covars2 = c("age","BMI","PMI","Smoking")
names(covars2) <- covars2
checkList2 <- lapply(c("pc3 only",covars2), function(x) filter(la.all, adj.P.Val<.05, covar==x)$gene_id)
names(checkList2) <- c("pc3", covars2)


## t test
t1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                     results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(covar="none")

t2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                     results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(covar="pc3 only")

t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                     results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))

la.t.all = bind_rows(t1, t2, t3) %>% 
  mutate(coef= factor(coef, levels=comparisons),
         covar=factor(covar, levels=covars))

### compare number of DEGs
degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(la.all, covar==y, adj.P.Val<.05)$gene_id
  filter(la.t.all, covar==y, adj.P.Val<.05, gene_id %in% sig.genes)
})) 

p2 <- ggplot(group_by(degs, covar, coef, .drop=F) %>% tally(), 
         aes(x=coef, y=n, color=covar))+
  ggbeeswarm::geom_beeswarm(size=2)+scale_color_manual(values=covars.colors)+
  scale_x_discrete(labels=gsub("_", "\n", comparisons))+
  labs(title="Layer-adjusted DEGs (t test <.05 and F test <.05)", subtitle=results_set,
       y="# genes", x="")+
  theme_minimal()


### volcano with F test filter
xlims = ceiling(max(abs(la.t.all$logFC)))
if(max(abs(la.t.all$logFC))<(xlims-.5)) xlims=xlims-.5

v1 <- ggplot(filter(la.t.all, covar %in% covars1), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3, color="grey50")+
  geom_point(data=filter(degs, covar %in% covars1), size=.3, color="red2")+
  facet_grid(rows=vars(coef), cols=vars(covar))+
  geom_vline(aes(xintercept=-.3), lty=2)+geom_vline(aes(xintercept=.3), lty=2)+
  xlim(-xlims,xlims)+
  labs(title= paste0("Layer-adjusted t test (", results_set, ")"),
       subtitle="(red DEGs are F test filtered)")+
  theme_bw()+theme(panel.grid.minor=element_blank())

v2 <- ggplot(filter(la.t.all, covar %in% covars2), aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point(size=.3, color="grey50")+
  geom_point(data=filter(degs, covar %in% covars2), size=.3, color="red2")+
  facet_grid(rows=vars(coef), cols=vars(covar))+
  geom_vline(aes(xintercept=-.3), lty=2)+geom_vline(aes(xintercept=.3), lty=2)+
  xlim(-xlims,xlims)+
  labs(title= paste0("Layer-adjusted t test (", results_set, ")"),
       subtitle="(red DEGs are F test filtered)")+
  theme_bw()+theme(panel.grid.minor=element_blank())

pdf(file=paste0("plots/07-1_covariate_sensitivity/layer-adjusted_", results_set, 
                "_covariate-comparisons.pdf"))
p1
upset(fromList(checkList1), text.scale=1.5)
grid.text("QC covariates", x = 0.65, y = 0.95, gp = gpar(fontsize = 16))
grid.text("L-A F test adj p<.05", x = 0.65, y = 0.9, gp = gpar(fontsize = 12))
upset(fromList(checkList2), text.scale=1.5)
grid.text("Donor covariates", x = 0.65, y = 0.95, gp = gpar(fontsize = 16))
grid.text("L-A F test adj p<.05", x = 0.65, y = 0.9, gp = gpar(fontsize = 12))
p2
grid.arrange(rasterize(v1, dpi=250), top="QC covariates")
grid.arrange(rasterize(v2, dpi=250), top="Donor covariates")
dev.off()

cat("\n\nL-A plots saved to:",paste0("plots/07-1_covariate_sensitivity/layer-adjusted_", results_set,
                "_covariate-comparisons.pdf"),"\n\n")


# layer-restricted
## F test
f1.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                         results_set, "_rev-gene-input_F-test_all.csv")) %>%
  mutate(covar="none")

f2.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                         results_set, "_rev-gene-input_F-test_all.csv")) %>%
  mutate(covar="pc3 only")

f3.all = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                         results_set, "_rev-gene-input_F-test_all_all-covar-models.csv"))

lr.all = bind_rows(f1.all, f2.all, f3.all) %>%
  mutate(covar=factor(covar, levels=covars))

#filter(lr.all, adj.P.Val<.05) %>% group_by(covar) %>% tally()
#filter(lr.all, adj.P.Val<.05) %>% group_by(covar) %>% slice_min(F) %>%
#  select(covar, F, P.Value, adj.P.Val)

p1 <- ggplot(filter(lr.all, adj.P.Val<.05), aes(x=covar, fill=covar))+
  geom_bar(stat="count", color="black", linewidth=.5)+
  geom_hline(aes(yintercept=nrow(filter(lr.all, adj.P.Val<.05, covar=="pc3 only"))), lty=2)+
  scale_fill_manual(values=covars.colors, guide="none")+
  labs(y="# genes", title="Layer-restricted F test adj p<.05", subtitle=results_set,
       x="covariate in model")+
  theme_minimal()


checkList1 <- lapply(c("pc3 only", covars1), function(x) filter(lr.all, adj.P.Val<.05, covar==x)$gene_id)
names(checkList1) <- c("pc3", covars1)

checkList2 <- lapply(c("pc3 only",covars2), function(x) filter(lr.all, adj.P.Val<.05, covar==x)$gene_id)
names(checkList2) <- c("pc3", covars2)


## t test
t1 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                     results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(covar="none")

t2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                     results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
  mutate(covar="pc3 only")

t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                     results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))

lr.t.all = bind_rows(t1, t2, t3) %>% 
  mutate(sex.group= factor(paste(sex, group, sep="_"), levels=comparisons),
         covar=factor(covar, levels=covars),
         cluster=factor(cluster, levels=names(comp_names)))

### compare number of DEGs
degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(lr.all, covar==y, adj.P.Val<.05)$gene_id
  filter(lr.t.all, covar==y, adj.P.Val<.05, gene_id %in% sig.genes)
}))

p2 <-ggplot(group_by(degs, covar, sex.group, cluster, .drop=F) %>% tally(), 
       aes(x=sex.group, y=n, color=covar))+
  ggbeeswarm::geom_beeswarm(size=2)+scale_color_manual(values=covars.colors)+
  facet_wrap(vars(cluster), scales="free")+
  labs(title="Layer-restricted DEGs (t test <.05 and F test <.05)", 
       subtitle=results_set, y="# genes", x="")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        panel.border = element_rect(fill=NA, color="grey70"),
                        strip.background = element_rect(fill="grey80", color=NA))


### volcano with F test filter
vlist1 <- lapply(comparisons, function(x) {
  sg.df = filter(lr.t.all, covar %in% covars1, sex.group==x)
  xlims = ceiling(max(abs(sg.df$logFC)))
  if(max(abs(sg.df$logFC))<(xlims-.5)) xlims=xlims-.5
  tmp = do.call(rbind, lapply(covars1, function(y) {
    sig.genes = filter(f3.all, covar==y, adj.P.Val<.05)$gene_id
    filter(sg.df, covar==y, adj.P.Val<.05, gene_id %in% sig.genes)
  }))
  pp <- ggplot(sg.df, aes(x=logFC, y=-log10(adj.P.Val)))+
    geom_point(size=.3, color="grey50")+geom_point(data=tmp, size=.3, color="red2")+
    facet_grid(rows=vars(cluster), cols=vars(covar))+
    geom_vline(aes(xintercept=-.3), lty=2)+geom_vline(aes(xintercept=.3), lty=2)+
    xlim(-xlims, xlims)+
    labs(title= paste0("Layer-restricted t test (", results_set, "): QC covariates"),
         subtitle="(red DEGs are F test filtered)")+
    theme_bw()+theme(panel.grid.minor=element_blank())
  rasterize(pp, dpi=250)
})

vlist2 <- lapply(comparisons, function(x) {
  sg.df = filter(lr.t.all, covar %in% covars2, sex.group==x)
  xlims = ceiling(max(abs(sg.df$logFC)))
  if(max(abs(sg.df$logFC))<(xlims-.5)) xlims=xlims-.5
  tmp = do.call(rbind, lapply(covars2, function(y) {
    sig.genes = filter(f3.all, covar==y, adj.P.Val<.05)$gene_id
    filter(lr.t.all, covar==y, sex.group==x, adj.P.Val<.05, gene_id %in% sig.genes)
  }))
  pp <- ggplot(sg.df, aes(x=logFC, y=-log10(adj.P.Val)))+
    geom_point(size=.3, color="grey50")+geom_point(data=tmp, size=.3, color="red2")+
    facet_grid(rows=vars(cluster), cols=vars(covar))+
    geom_vline(aes(xintercept=-.3), lty=2)+geom_vline(aes(xintercept=.3), lty=2)+
    xlim(-xlims, xlims)+
    labs(title= paste0("Layer-restricted t test (", results_set, "): Donor covariates"),
         subtitle="(red DEGs are F test filtered)")+
    theme_bw()+theme(panel.grid.minor=element_blank())
  rasterize(pp, dpi=250)
})

test = c(vlist1, vlist2)
test = test[c(1,7,2,8,3,9,4,10,5,11,6,12)]

pdf(file=paste0("plots/07-1_covariate_sensitivity/layer-restricted_", results_set, 
                "_covariate-comparisons.pdf"))
p1
upset(fromList(checkList1), text.scale=1.5)
grid.text("QC covariates", x = 0.65, y = 0.95, gp = gpar(fontsize = 16))
grid.text("L-R F test adj p<.05", x = 0.65, y = 0.9, gp = gpar(fontsize = 12))
upset(fromList(checkList2), text.scale=1.5)
grid.text("Donor covariates", x = 0.65, y = 0.95, gp = gpar(fontsize = 16))
grid.text("L-R F test adj p<.05", x = 0.65, y = 0.9, gp = gpar(fontsize = 12))
p2
marrangeGrob(test, ncol=1, nrow=1, top = quote(rep(comparisons, each=2)[[g]]))
dev.off()

cat("\n\nL-R plots saved to:", paste0("plots/07-1_covariate_sensitivity/layer-restricted_", results_set,
                "_covariate-comparisons.pdf"),"\n\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
