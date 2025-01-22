# Sample info sources
### SOURCE 1
"DLPFC cross-disorders brain collection bookkeeping updated 1/15/2025" is a google sheet file (URL: https://docs.google.com/spreadsheets/d/1hBW74HWz-BqHbPEROHsSOZhRYizHDTSPY9l3_RcerSY/edit?gid=0#gid=0). Values from rows 4 through 225 and columns A (brain number), M (diagnosis), N (age), O (sex), P (PMI), Q (RIN), and Z (Study 2 MDD/BPD) from the tab named "Big_240_DLPFC_dissections" were copied to an offline spreadsheet. These data were then filtered only to samples used in the 'Study 2 MDD/BPD', the brain number column was cleaned to separate tissue notes from the brain number, and diagnoses were relabeled to NTC, MDD, BPD. These data were copied into a clean csv that was transferred to JHPCE and named "DLPFC_cross-disorders_demographics_MBv.csv"

### SOURCE 2
"Visium_DATA_2025-01-22_1406.csv" was pulled from REDCap on January 22, 2025. This data sheet contains the donor brain number, Visium slide codes (e.g., V13B23-XXX), and array positions (e.g., A1) for many samples including the MBv samples.

# Sequencing rounds
Round 1: An initial sequencing run was performed on the first 6 non-pilot slides (V13F27-338, V13F27-348, V13Y10-020, V13Y10-021, V13Y10-022, V13Y10-023). 
Round 2: The remaining 34 slides were run together (round 2). An experimental error occurred on V13B23-339_A1 that was caught prior to submitting the samples for sequencing. This sample was not sequenced and a new tissue section and sample preparation was performed on V13B23-283_B1 which was submitted with other Visium slides for a lateral septum experiment.
Round 3: Quality control performed on the sequenced data from all 40 slides revealed two slides with substantially fewer reads across all samples (one from round 1 V13Y10-020 and one from round 2 V13B23-331). See "plots/02_build_spe/r1_r2_spaceranger_overview.pdf". These new tissue sections were selected for the 8 samples on these slides. V13Y10-020 --> V13B23-282, V13B23-331 --> V13B23-279. 

# Compiling demographic information for MBv project
See script found in "code/02_build_spe/01_raw_spe_hdf5.r" for how the sample info sources were combined. Following the details listed above, specific filters were applied when compiling the sample info and generating the file list used as input for constructing the spe object: 
-- removal of slides V13Y10-020 and V13B23-331 from input list and metadata data frame (see sequencing round 3 notes)
-- removal of sample V13B23-339_A1 from metadata data frame (sequencing round 2 notes)
Lastly, there is one brain donor (Br5366) which provided sections for two different samples in this experiment: V13F27-338_C1 (round 1) and V13B23-380_A1 (round 2). This means that although there are 120 samples, there are only 119 donors and that the MDD male group only has 19 donors while all of the other diagnosis*sex groups have 20. 

