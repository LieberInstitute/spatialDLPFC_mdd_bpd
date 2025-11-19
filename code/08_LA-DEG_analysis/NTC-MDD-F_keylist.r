keylist <- list()
#reactome smoothed
keylist[["smoothed"]][["Reactome"]] = 
  list("Ribo"=c("Eukaryotic Translation Elongation R-HSA-156842",
                "Nonsense Mediated Decay (NMD) Enhanced By Exon Junction Complex (EJC) R-HSA-975957",
                "Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency R-HSA-9633012",
                "Cellular Response To Starvation R-HSA-9711097"),
       "Inter"=c("Cytokine Signaling In Immune System R-HSA-1280215"),
       "Complement"=c("Complement Cascade R-HSA-166658"),
       "Platelet"=c("Hemostasis R-HSA-109582","Platelet Degranulation R-HSA-114608"),
       "Neuro"=c("GABA Synthesis, Release, Reuptake And Degradation R-HSA-888590",
                 "Neuronal System R-HSA-112316")
  )

#go-bp smoothed
keylist[["smoothed"]][["GO-BP"]] = 
  list("Copper"=c("Negative Regulation of Growth (GO:0045926)"),
       "Cytokine"=c("Regulation of Cytokine Production (GO:0001817)",
                    "Cytokine-Mediated Signaling Pathway (GO:0019221)",
                    "Response to Cytokine (GO:0034097)",
                    "Cellular Response to Cytokine Stimulus (GO:0071345)",
                    "Regulation of Canonical NF-kappaB Signal Transduction (GO:0043122)"),
       "Inter"=c(),
       "Apop"=c("Negative Regulation of Apoptotic Process (GO:0043066)",
                "Positive Regulation of Apoptotic Process (GO:0043065)"),
       "SMAD"=c("SMAD Protein Signal Transduction (GO:0060395)"),
       "Translation"=c("Cytoplasmic Translation (GO:0002181)"),
       "Phago"=c("Regulation of Phagocytosis (GO:0050764)"),
       "Inflamm"=c("Inflammatory Response (GO:0006954)"),
       "mRNA"=c("Positive Regulation of mRNA Catabolic Process (GO:0061014)",
                "3'-UTR-mediated mRNA Stabilization (GO:0070935)",
                "Nuclear-Transcribed mRNA Catabolic Process, Nonsense-Mediated Decay (GO:0000184)"),
       "Synap"=c("Chemical Synaptic Transmission (GO:0007268)",
                 "Modulation of Chemical Synaptic Transmission (GO:0050804)",
                 "Synapse Assembly (GO:0007416)")
  )

#go-cc smoothed
keylist[["smoothed"]][["GO-CC"]] = 
  list("Adhesion"=c("Focal Adhesion (GO:0005925)"),
       "Ribo"=c("Large Ribosomal Subunit (GO:0015934)",
                "Small Ribosomal Subunit (GO:0015935)"),
       "Platelet"=c("Platelet Alpha Granule (GO:0031091)"),
       "Secretory"=c("Secretory Granule Membrane (GO:0030667)", 
                     "Secretory Granule Lumen (GO:0034774)",
                     "Cytoplasmic Vesicle Lumen (GO:0060205)"),
       "Ficolin"=c("Ficolin-1-Rich Granule (GO:0101002)")
  )

#reactome seurat
keylist[["seurat"]][["Reactome"]] = 
  list("Ribo"=c("Eukaryotic Translation Elongation R-HSA-156842",
                "Nonsense Mediated Decay (NMD) Enhanced By Exon Junction Complex (EJC) R-HSA-975957",
                "Response Of EIF2AK4 (GCN2) To Amino Acid Deficiency R-HSA-9633012",
                "Cellular Response To Starvation R-HSA-9711097"),
       #"Cytokine"=c("Cytokine Signaling In Immune System R-HSA-1280215"),
       "Inter"=c("Cytokine Signaling In Immune System R-HSA-1280215"),
       "Complement"=c("Complement Cascade R-HSA-166658"),
       "Platelet"=c("Hemostasis R-HSA-109582","Platelet Degranulation R-HSA-114608"),
       "Neuro"=c("GABA Synthesis, Release, Reuptake And Degradation R-HSA-888590",
                 "Neuronal System R-HSA-112316"),
       "TLR"=c("Toll-like Receptor Cascades R-HSA-168898",
               "ER-Phagosome Pathway R-HSA-1236974")
  )

#go-bp seurat
keylist[["seurat"]][["GO-BP"]] = 
  list("Copper"=c("Negative Regulation of Growth (GO:0045926)",
                  "Intracellular Monoatomic Cation Homeostasis (GO:0030003)",
                  "Cellular Response to Metal Ion (GO:0071248)"),
       "Cytokine"=c("Regulation of Cytokine Production (GO:0001817)",
                    "Cytokine-Mediated Signaling Pathway (GO:0019221)",
                    "Response to Cytokine (GO:0034097)",
                    "Cellular Response to Cytokine Stimulus (GO:0071345)"),
       "Inter"=c(),
       "Apop"=c("Negative Regulation of Apoptotic Process (GO:0043066)",
                "Positive Regulation of Apoptotic Process (GO:0043065)"),
       "Neutrophil"=c("Granulocyte Chemotaxis (GO:0071621)"),
       "Phago"=c("Regulation of Phagocytosis (GO:0050764)"),
       "Inflamm"=c("Inflammatory Response (GO:0006954)",
                   "Positive Regulation of Inflammatory Response (GO:0050729)"),
       "Virus|Viral"=c("Defense Response to Virus (GO:0051607)",
                       "Negative Regulation of Viral Process (GO:0048525)"),
       "Bacterium"=c(),
       "mRNA"=c("Positive Regulation of mRNA Catabolic Process (GO:0061014)",
                "Nuclear-Transcribed mRNA Catabolic Process (GO:0000956)"),
       "Blood"=c(),
       "Lipo"=c(),
       "Transcription Factor"=c("Positive Regulation of DNA-binding Transcription Factor Activity (GO:0051091)"),
       "Synap"=c("Chemical Synaptic Transmission (GO:0007268)",
                 "Modulation of Chemical Synaptic Transmission (GO:0050804)")
  )

#go-cc seurat
keylist[["seurat"]][["GO-CC"]] = 
  list("Adhesion"=c("Focal Adhesion (GO:0005925)"),
       "Ribo"=c("Large Ribosomal Subunit (GO:0015934)",
                "Small Ribosomal Subunit (GO:0015935)"),
       "Platelet"=c("Platelet Alpha Granule (GO:0031091)"),
       "Secretory"=c("Secretory Granule Membrane (GO:0030667)", 
                     "Secretory Granule Lumen (GO:0034774)",
                     "Cytoplasmic Vesicle Lumen (GO:0060205)"),
       "Ficolin"=c("Ficolin-1-Rich Granule (GO:0101002)")
  )
