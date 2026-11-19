# About
Principal demographic data is derived from two major sources (see below). Throughout the analysis process,
additional donor information was explored and is present in smaller distinct CSV files. 

<br>
>[!WARNING]
>This data is for internal LIBD reference only and is not suitable for public distribution due to the 
>abundance of sensitive PHI that could be used to identify donors. Donor demographic info approved
>for public distribution references these files. The build info for the public files can be found in
>the corresponding `code` files (e.g., Globus, publication Tables and Supplemental Tables).
<br>

In addition to metadata CSV files, this directory contains two R scripts. The code in `getMBvSampleInfo_function.R` 
is a wrapper to build the core demographic info from the principal sources. This function is 
utilized in the initial SPE build instructions and in the `build_demographics_table.R` script to
generate the `n119_demographics.csv` present in this directory. The `n119_demographics.csv` will 
serve as a reference for all of the scripts constructing any tables or figure/SF panels present in
the manuscript that require the use of PHI.

# Sample info sources
### Principal Source 1
"DLPFC cross-disorders brain collection bookkeeping updated 1/15/2025" is a google sheet file 
(URL: https://docs.google.com/spreadsheets/d/1hBW74HWz-BqHbPEROHsSOZhRYizHDTSPY9l3_RcerSY/edit?gid=0#gid=0). 
Values from rows 4 through 225 and columns A (brain number), M (diagnosis), N (age), O (sex), P (PMI), Q (RIN), 
and Z (Study 2 MDD/BPD) from the tab named "Big_240_DLPFC_dissections" were copied to an offline spreadsheet. 
These data were then filtered only to samples used in the 'Study 2 MDD/BPD', the brain number column was cleaned 
to separate tissue notes from the brain number, and diagnoses were relabeled to NTC, MDD, BPD. These data were 
copied into a clean csv that was transferred to JHPCE and named "DLPFC_cross-disorders_demographics_MBv.csv"

### Principal Source 2
"Visium_DATA_2025-01-22_1406.csv" was pulled from REDCap on January 22, 2025. This data sheet contains the donor
brain number, Visium slide codes (e.g., V13B23-XXX), and array positions (e.g., A1) for many samples including 
the MBv samples.
UPDATE:  During sample prep samples on slide 301, the samples were mounted in reverse intended order (instead 
of A-D they were mounted D-A). In May 2025, sex chr sanity checks revealed that the actual sample order was not 
reflected in REDCap. SCP corrected REDCap and JRT modified `getMBvSampleInfo_function.R` to correct the metadata.

### Additional Sources
NeuroPath provided several additional files in response to a couple different requests during the course of
analysis. 
- Page_BMI_Smoking_CrossDisorder_120825.csv: BMI, next-of-kin reported smoking status, nicotine metabolite
measurements (if collected)
- Bipolar_subtype.csv: not formally included in analysis but used to spot-check seemingly bimodal BD DEGs
(concluded not related to BD subtype) and used to check for PRS differences calculated from GWAS-BD and GWAS-SCZ
- MBv_suicide.csv: Group summary statistics for suicide COD

# Additional Notes
### Sequencing rounds
Round 1: An initial sequencing run was performed on the first 6 non-pilot slides (V13F27-338, V13F27-348, V13Y10-020, 
V13Y10-021, V13Y10-022, V13Y10-023). These were actually run for different cores (see batches).
Round 2: The remaining 34 slides were run together (round 2). An experimental error occurred on V13B23-339_A1 that 
was caught prior to submitting the samples for sequencing. This sample was not sequenced and a new tissue section 
and sample preparation was performed on V13B23-283_B1 which was submitted with other Visium slides for a lateral 
septum experiment.
Round 3: Quality control performed on the sequenced data from all 40 slides revealed two slides with substantially 
fewer reads across all samples (one from round 1 V13Y10-020 and one from round 2 V13B23-331). See 
"plots/01_spaceranger/INTERNAL_r1-r2-r3_spaceranger-overview.png". These new tissue sections were selected for 
the 8 samples on these slides. V13Y10-020 --> V13B23-282, V13B23-331 --> V13B23-279. 

### Repeated MDD male Compiling demographic information for MBv project
One brain donor (Br5366) provided sections for two different samples in this experiment: V13F27-338_C1 (round 1) 
and V13B23-380_A1 (round 2). This means that although there are 120 samples, there are only 119 donors and that 
the MDD male group only has 19 donors while all of the other diagnosis*sex groups have 20. One of these samples
(V13F27-338_C1) was dropped after QC.

