#This will be the app file for this project
library("spatialLIBD")
library("here")
library("devtools")
library("HDF5Array")
library("markdown")
library("SingleCellExperiment")
library("iSEE")
## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())


## Load the data
posit_connect_file <- "/r_data/lcollado/Posit_Connect_shiny_apps/dlpfc/DLPFC_MBv_pseudobulk/iSEE_pseudobulk-spe_both-annotations.rds"


if (file.exists(posit_connect_file)) {
    ## Location for the https://conn1.libd.org/ server
    load(posit_connect_file)
} else {
    spe <- readRDS(here::here("dlpfc", "DLPFC_MBv_pseudobulk", "processed-data", "iSEE_pseudobulk-spe_both-annotations.rds"))
}

# replace rownames
rownames(spe) <- rowData(spe)$gene_name
rownames(metadata(spe)[[1]]) <- metadata(spe)[[1]]$gene_name
rownames(metadata(spe)[[2]]) <-	metadata(spe)[[2]]$gene_name

# modify coldata factors (in future can just re-save object with desired factors)
spe$dx.sex = factor(paste(spe$sex, spe$condition), levels=c("F NTC","F MDD","F BPD","M NTC","M MDD","M BPD"),
                    labels=c("F NTC","F MDD","F BD","M NTC","M MDD","M BD"))
spe$diagnosis = factor(spe$condition, levels=c("NTC","MDD","BPD"), labels=c("NTC","MDD","BD"))
spe$domain_annotation = factor(paste0(spe$annotation, ": ", spe$domain),
                           levels=c(paste0("domain-SP: ", c("L1","L2","L3.4","L5","L6","WM")), 
                                    paste0("domain-CT: ", c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))),
                           labels=c("domain-SP: L1", "domain-SP: L2", "domain-SP: L3/4", "domain-SP: L5", "domain-SP: L6",  "domain-SP: WM",
                                    "domain-CT: M/V", "domain-CT: Ast", "domain-CT: L2/3", "domain-CT: L4", "domain-CT: Inb", "domain-CT: L5", "domain-CT: L6", "domain-CT: Olg"))
spe$consensus_domains = factor(spe$domain, levels=c("Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4","Inhb","L5","L6","WM","Oligo"),
                              labels=c("M/V","Ast","L1","L2/L3","L2/L3","L3/L4","L3/L4","Inb","L5","L6","WM/O","WM/O"))
spe$annotation = factor(spe$annotation, levels=c("domain-SP","domain-CT"))

# to make it easiest I might just filter spe to genes run in both models
spe2 <- spe[intersect(rownames(spe@metadata[[1]]), rownames(spe@metadata[[2]])),]
spe2@metadata[[1]] = spe2@metadata[[1]][rownames(spe2),]
spe2@metadata[[2]] = spe2@metadata[[2]][rownames(spe2),]

# subset to only necessary colData so less confusing plotting options for users
colData(spe2) <- colData(spe2)[,c("sex","diagnosis","dx.sex","annotation","domain_annotation","consensus_domains")]

# read in initial and custom functions
source("custom_table_functions.R")
source("initial.R", print.eval = TRUE)

# change panel box color
spe2 <- registerAppOptions(spe2, panel.color=c(RowDataTable="#000000",
                                               FeatureAssayPlot="#99999C",
                                               DomainRestrictedDE="#99999C",
                                               WholeTissueDE="#99999C"))

rowData(spe2) <- addDEGsToRowData(spe2)

iSEE(spe2, initial=initial, 
     colormap=ecm,
     appTitle = "pseudobulk dlPFC spatial data: NTC/MDD/BD",
     tour=mytour)
