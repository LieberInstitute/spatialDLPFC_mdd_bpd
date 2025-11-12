keylist <- list()
#reactome smoothed
keylist[["smoothed"]][["Reactome"]] = 
  list("Collagen"=c("Collagen Formation R-HSA-1474290","Extracellular Matrix Organization R-HSA-1474244"),
       "Respiratory"=c("Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517"),
       "GPCR"=c("GPCR Ligand Binding R-HSA-500792"),
       "NTRK"=c("Signaling By NTRK1 (TRKA) R-HSA-187037"),
       "RNA"=c("Eukaryotic Translation Elongation R-HSA-156842",
               "Nonsense Mediated Decay (NMD) Enhanced By Exon Junction Complex (EJC) R-HSA-975957",
               "rRNA Processing R-HSA-72312","Influenza Infection R-HSA-168255"))
#go-bp smoothed
keylist[["smoothed"]][["GO-BP"]] = 
  list("Proton"=c("Cellular Respiration (GO:0045333)","Proton Motive Force-Driven ATP Synthesis (GO:0015986)"),
       "Extracellular"=c("Extracellular Matrix Organization (GO:0030198)"),
       "Translation"=c("Translation (GO:0006412)","Gene Expression (GO:0010467)"))

#gocc smoothed
keylist[["smoothed"]][["GO-CC"]] = 
  list("Ribosomal"=c("Large Ribosomal Subunit (GO:0015934)","Small Ribosomal Subunit (GO:0015935)"))


#reactome seurat
keylist[["seurat"]][["Reactome"]] = 
  list("Respiratory"=c("Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517"),
       "GPCR"=c("GPCR Ligand Binding R-HSA-500792"),
       "NTRK"=c("Signaling By NTRK1 (TRKA) R-HSA-187037"),
       "RNA"=c("Eukaryotic Translation Elongation R-HSA-156842",
               "Nonsense Mediated Decay (NMD) Enhanced By Exon Junction Complex (EJC) R-HSA-975957",
               "rRNA Processing R-HSA-72312","Influenza Infection R-HSA-168255"))

#go-bp seurat
keylist[["seurat"]][["GO-BP"]] = 
  list("Copper"=c("Response to Copper Ion (GO:0046688)"),
       "Proton"=c("Cellular Respiration (GO:0045333)",
                  "Proton Motive Force-Driven ATP Synthesis (GO:0015986)"))

#go-cc seurat
keylist[["seurat"]][["GO-CC"]] = 
  list("Ribosomal"=c("Large Ribosomal Subunit (GO:0015934)","Small Ribosomal Subunit (GO:0015935)"),
       "Vesicle"=c("Extracellular Vesicle (GO:1903561)","Clathrin-Coated Vesicle Membrane (GO:0030665)"))
