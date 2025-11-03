
keylist <- list()

#seurat
#Reactome
keylist[["seurat"]][["Reactome"]] = 
  list("Electron"=c("Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517"),
       "Mitochondrial"=c("Mitochondrial Translation R-HSA-5368287","Translation R-HSA-72766"),
       "GPCR"=c("GPCR Ligand Binding R-HSA-500792","GPCR Downstream Signaling R-HSA-388396"),
       "Complement"=c("Complement Cascade R-HSA-166658","Cell Recruitment (Pro-Inflammatory Response) R-HSA-9664424"),
       "Platelet"=c("Platelet Degranulation R-HSA-114608"),
       "Interleukin"=c("Cytokine Signaling In Immune System R-HSA-1280215")
  )


#GO-BP
keylist[["seurat"]][["GO-BP"]] = 
  list("Mitochondrial"=c("Cellular Respiration (GO:0045333)","Mitochondrial Respiratory Chain Complex Assembly (GO:0033108)",
                         "Proton Motive Force-Driven ATP Synthesis (GO:0015986)",
                         "Mitochondrial Translation (GO:0032543)","Translation (GO:0006412)"),
       "Micro"=c("Microglial Cell Activation (GO:0001774)","Response to Axon Injury (GO:0048678)"),
       "Chemo"=c("Positive Regulation of Chemotaxis (GO:0050921)","Regulation of Chemokine Production (GO:0032642)"),
       "Leukocyte"=c("Leukocyte Adhesion to Vascular Endothelial Cell (GO:0061756)"),
       "Antigen"=c("Peptide Antigen Assembly With MHC Protein Complex (GO:0002501)",
                   "Positive Regulation of Immune Response (GO:0050778)"),
       "Cytokine"=c("Response to Cytokine (GO:0034097)","Cellular Response to Cytokine Stimulus (GO:0071345)"),
       "Inflamm"=c("Regulation of Inflammatory Response (GO:0050727)"),
       "Localization"=c("Positive Regulation of Protein Localization to Nucleus (GO:1900182)"),
       "Stabilization"=c("mRNA Stabilization (GO:0048255)"),
       "Toll"=c(),
       "Copper"=c("Copper Ion Homeostasis (GO:0055070)")
  )

#GO-CC
keylist[["seurat"]][["GO-CC"]] = 
  list("Vesicle Membrane"=c("Clathrin-Coated Vesicle Membrane (GO:0030665)","MHC Protein Complex (GO:0042611)",
                            "ER to Golgi Transport Vesicle Membrane (GO:0012507)"),
       "Platelet"=c("Platelet Alpha Granule (GO:0031091)","Secretory Granule Lumen (GO:0034774)")
  )

#smoothed
#Reactome
keylist[["smoothed"]][["Reactome"]] = 
  list("Electron"=c("Respiratory Electron Transport R-HSA-611105"),
       "Mitochondrial"=c("Mitochondrial Translation R-HSA-5368287"),
       "Interleukin"=c("Cytokine Signaling In Immune System R-HSA-1280215"),
       "Complement"=c("Complement Cascade R-HSA-166658"),
       "Alpha"= c("G Alpha (I) Signaling Events R-HSA-418594"),
       "FOX"=c("FOXO-mediated Transcription Of Cell Cycle Genes R-HSA-9617828"),
       "Platelet"=c("Platelet Degranulation R-HSA-114608")
  )

#GO-BP
keylist[["smoothed"]][["GO-BP"]] = 
  list("Mitochondrial"=c("Cellular Respiration (GO:0045333)","Mitochondrial Respiratory Chain Complex Assembly (GO:0033108)",
                         "Mitochondrial Translation (GO:0032543)"),
       #"Interleukin"=c("Microglial Cell Activation (GO:0001774)","Phagocytosis (GO:0006909)",
       #                 "Cellular Response to Tumor Necrosis Factor (GO:0071356)"),
       "Chemo"=c("Microglial Cell Activation (GO:0001774)","Phagocytosis (GO:0006909)",
                 "+ Reg of Phosphatidylinositol 3-Kinase/Prot Kinase B Signal Transduction (GO:0051897)"),
       "Cellular Response"=c("Cellular Response to Transforming Growth Factor Beta Stimulus (GO:0071560)" ,
                             "Cellular Response to Tumor Necrosis Factor (GO:0071356)"),
       "Proliferation"=c("Positive Regulation of Immune Response (GO:0050778)"),
       "Differentiation"=c("Myeloid Leukocyte Differentiation (GO:0002573)","Regulation of Epithelial Cell Differentiation (GO:0030856)"),
       "Leukocyte"=c("Leukocyte Adhesion to Vascular Endothelial Cell (GO:0061756)"),
       "Antigen"=c("Peptide Antigen Assembly With MHC Protein Complex (GO:0002501)"),
       "Stabilization"=c("mRNA Stabilization (GO:0048255)"),
       "STAT"=c("Regulation of Tyrosine Phosphorylation of STAT Protein (GO:0042509)"),
       "Development"=c("Central Nervous System Development (GO:0007417)")
  )


#GO-CC
keylist[["smoothed"]][["GO-CC"]] = 
  list("Vesicle Membrane"=c("Clathrin-Coated Vesicle Membrane (GO:0030665)","MHC Protein Complex (GO:0042611)",
                            "ER to Golgi Transport Vesicle Membrane (GO:0012507)"),
       "Granule"=c("Platelet Alpha Granule (GO:0031091)","Secretory Granule Membrane (GO:0030667)")
  )

