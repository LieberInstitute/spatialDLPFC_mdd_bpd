#plot ORA results
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

results_set = "smoothed-k9-1663"
#results_set = "seurat-pc30"

ora <- readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-adjusted_", 
                      results_set,"_F-test-padj05_Reactome-ORA.rda"))
#smoothed
la.lr = c("MAPK Targets/ Nuclear Events Mediated By MAP Kinases",
          "Cellular Response To Starvation",
          "Peptide Ligand-Binding Receptors",
          "Cytokine Signaling In Immune System",
          "Chaperone Mediated Autophagy"
          )
la.only = c("Respiratory Electron Transport",
            "GABA Synthesis, Release, Reuptake And Degradation",
            "Complement Cascade",
            "Axon Guidance")
lr.only = c("Signaling By NTRKs",
            "Toll-like Receptor Cascades",
            "Oxidative Stress Induced Senescence",
            "FOXO-mediated Transcription Of Cell Cycle Genes")
#seurat
la.lr = c("MAPK Targets/ Nuclear Events Mediated By MAP Kinases",
          "Signaling By NTRKs",
          "Peptide Ligand-Binding Receptors",
          "Cytokine Signaling In Immune System",
          "Chaperone Mediated Autophagy")
la.only = c("Cellular Response To Starvation",
            "Mitochondrial Biogenesis",
            "Respiratory Electron Transport",
            "TP53 Regulates Metabolic Genes",
            "Complement Cascade",
            "GABA Synthesis, Release, Reuptake And Degradation",
            "Axon Guidance")
lr.only = c("Toll-like Receptor Cascades",
            #"Spry Regulation Of FGF Signaling",
            "ER-Phagosome Pathway",
            "Suppression Of Phagosomal Maturation",
            "Circadian Clock",
            "Cellular Senescence",
            "FOXO-mediated Transcription Of Cell Cycle Genes")

filter(ora@result, Description %in% la.lr)
filter(ora@result, Description %in% la.only)
filter(ora@result, Description %in% lr.only)

ora_la = ora
ora_la@result <- ora_la@result[ora_la@result$Description %in% c(la.lr, la.only),]
dotplot(ora_la, x="Count", showCategory=nrow(ora_la@result))+
  scale_size_continuous(range=c(1,6), limits=c(0,.2), breaks=c(0,.1,.2))+
  xlim(0,45)+ggtitle("L-A (PRECAST)")

ora <- readRDS(paste0("processed-data/08_dx-sex_DEG_analysis/layer-restricted_", 
                      results_set,"_F-test-padj05_Reactome-ORA.rda"))

filter(ora@result, Description %in% la.lr)
filter(ora@result, Description %in% la.only)
filter(ora@result, Description %in% lr.only)

ora_lr = ora
ora_lr@result <- ora_lr@result[ora_lr@result$Description %in% c(la.lr, lr.only),]
dotplot(ora_lr, x="Count", showCategory=nrow(ora_lr@result))+
  scale_size_continuous(range=c(1,6), limits=c(0,.2), breaks=c(0,.1,.2))+
  xlim(0,30)+ggtitle("L-R (PRECAST)")
