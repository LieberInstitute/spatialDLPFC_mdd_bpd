setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("HDF5Array")
})

# Get sample info
source("code/01_spaceranger/INTERNAL_getMBvSampleInfo_function.r")
demo = getMBvSampleInfo(REDCapFile="Visium_DATA_2025-01-22_1406.csv",
                        demoFile="DLPFC_cross-disorders_demographics_MBv.csv")

# Identify samples with available FASTQ and spaceranger output
fastqs = grep("^V", list.files('raw-data/FASTQ'), value=T)
space = grep("^V", list.files("processed-data/01_spaceranger"), value=T)
## Check that the samples with FASTQ files and spaceranger output are identical
stopifnot(identical(fastqs, space))

# List any samples from demo data that are missing spaceranger output
missingSpaceranger = setdiff(demo$sample_id, space)
## One sample/donor was planned to run as V13B23-339_A1but had to be re-run as V13B23-283_B1
#### You can see Br8237 has two entries listed in demo table due to both intended and actual sample_ids being listed
cat("\nSamples from demo/exp data that are missing spaceranger output.\n(should only be V13B23-339_A1)\n")
missingSpaceranger
stopifnot(missingSpaceranger=="V13B23-339_A1")
demo = demo[demo$sample_id!="V13B23-339_A1",]
## Now that demo table and spaceranger output match, append outs path to demo table
stopifnot(identical(sort(space),sort(demo$sample_id)))
demo$sample_path = file.path("processed-data/01_spaceranger", demo$sample_id,"outs")

# Check for repeated samples for same donor
repeated_donors = demo$brnum[duplicated(demo$brnum)]
table(demo[demo$brnum %in% repeated_donors, c("slide")])
table(demo[demo$brnum %in% repeated_donors, c("seq","round")])
## Two slides were re-run due to lower quality, keep second run
cat("\nRemove samples from slides V13Y10-020 and V13B23-331 from the list used to generate the spe\n")
demo = demo[!demo$slide %in% c("V13Y10-020","V13B23-331"),]
## Now the only repeated donor is the MDD M accidentally included in the experiment twice
stopifnot(dim(demo)[1]==120)
cat("\nOnly repeated donor is known instance of double inclusion of same MDD male\n")
repeated_donors = demo$brnum[duplicated(demo$brnum)]
demo[demo$brnum %in% repeated_donors,]

#### In the initial build we moved forward with spe generation and ran QC for both instances
#### of donor Br5366. They were fairly equal. We decided to keep the sample that was run
#### in the large Psomagen batch and remove V13F27-338_C1. 
#### In the Nextflow/GEO build, only the Psomagen sample (V13B23-380_A1) is present for this donor.

