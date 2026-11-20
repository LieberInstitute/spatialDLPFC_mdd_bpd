args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
    stop("Usage: Rscript create_colorPalettes.R [output_rds]")
}
output_rds <- if (length(args) == 1L) args[[1L]] else "plots/colorPalettes.rds"

cpList <- list()

cpList$dx.pal = c("NTC"="#7f7f7f", "MDD"="#FB8861", "BD"="#9260b2") #sRGB friendly


cpList$low.res.light = c("Micro.Vasc"="#d05e46","Astro"="#F5D29E",
                         "L2"="#b3c7ac","L3"="#a3c7e4",
                         "L4"="#c2cfcf",
                         "L5"="#e5e8a0","L6"="#ef9e9f",
                         "Oligo"="#ede8e0","Inhb"="#b9a7c9")

cpList$low.res.bright =c("Micro.Vasc"="#911223","Astro"="#cfa45c",
                         "L2"="#5D9940", "L3"="#5095CD", 
                         "L4"="#c2cfcf", #"L4"="#85A0A0",
                         "L5"="#ddc94e","L6"="#E45C5F",
                         "Oligo"="#D1C4B0","Inhb"="#9377AC")

cpList$smoothed.light = c("L1"="#F5D29E", "L2"="#b3c7ac", 
                          "L3.4"="#a3c7e4", "L5"="#e5e8a0",
                          "L6"= "#ef9e9f", "WM"="#ede8e0")

cpList$smoothed.bright = c("L1"="#cfa45c", "L2"= "#5D9940",
                           "L3.4"= "#5095CD","L5"="#ddc94e",
                           "L6"="#E45C5F","WM"= "#D1C4B0")

cpList$transfer.light = c("M.V"="#d05e46","Ast"="#F5D29E",
                         "L2.3"="#83c7c7",
                         "L4"="#c2cfcf",
                         "L5"="#e5e8a0","L6"="#ef9e9f",
                         "Olg"="#ede8e0","Inb"="#b9a7c9")

cpList$transfer.bright =c("M.V"="#911223","Ast"="#cfa45c",
                         "L2.3"="#088F8F", 
                         "L4"="#c2cfcf", #"L4"="#85A0A0",
                         "L5"="#ddc94e","L6"="#E45C5F",
                         "Olg"="#D1C4B0","Inb"="#9377AC")

saveRDS(cpList, output_rds)
