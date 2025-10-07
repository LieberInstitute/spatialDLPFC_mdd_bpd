#look at plots produced by 03_evaluate-clusters and annotate the clusters
annotations <- list(
	precast_k5_1663 = c("Glia"=3, "L2"=2, "L3/4"=5, "L5"=1, "L6"=4),
	precast_k7_1663 = c("Vasc?"=4, "Vasc/L1"=3, "L2"=2, "L3/4"=6, "L5"=1, "L6"=5, "WM"=7),
	precast_k9_1663 = c("Vasc"=5, "L1"=3, "L2"=2, "L3/4"=7, "L5"=1, "L6"=6, "WM"=9, "GABA"=8, "low UMI"=4),
	precast_k12_1663 = c("Vasc"=12, "L1"=3, "L2 (1)"=2, "L2 (2)"=7, "L3/4 (1)"=9, "L3/4 (2)"=4,
		"L5"=1, "L6"=8, "WM (1)"=11, "WM (2)"=5, "GABA"=10, "low UMI"=6),
	precast_k7_HM = c("Vasc"=3, "L1"=1, "L2"=2, "L3/4"=5, "L5"=7, "L6"=4, "WM"=6),
	precast_k9_HM = c("Vasc"=3, "L1"=1, "L2"=2, "L3/4"=6, "L5"=9, "L6"=5, "WM"=8, "GABA"=7, "L2/3?"=4),
	precast_k9_1626 = c("Vasc"=5, "L1"=4, "L2"=9, "L3/4 (1)"=7, "L3/4 (2)"=2, "L5"=3, "L6"=1, "WM"=8, "low UMI"=6),
	precast_k7_1626 = c("L1"=4, "L2"=7, "L3/4"=2, "L5"=3, "L6"=1, "WM"=5, "low UMI"=6)
)

#based on earthy.pal
annot_colors <- list(
        precast_k5_1663 = c("Glia"="#911223", "L2"="#5D9940", 
		"L3/4"="#5095CD", "L5"="#ddc94e", "L6"="#E45C5F"),
        precast_k7_1663 = c("Vasc?"="#911223", "Vasc/L1"="#cfa45c", "L2"="#5D9940", 
		"L3/4"="#5095CD", "L5"="#ddc94e", "L6"="#E45C5F", "WM"="#D1C4B0"),
        precast_k9_1663 = c("Vasc"="#911223", "L1"="#cfa45c", "L2"="#5D9940", 
		"L3/4"="#5095CD", "L5"="#ddc94e", "L6"="#E45C5F", 
		"WM"="#D1C4B0", "GABA"="#9377AC", "low UMI"="grey"),
        precast_k12_1663 = c("Vasc"="#911223", "L1"="#cfa45c", 
		"L2 (1)"="#5D9940", "L2 (2)"="#b3c7ac", 
		"L3/4 (1)"="#5095CD", "L3/4 (2)"="#a3c7e4",
                "L5"="#ddc94e", "L6"="#E45C5F", 
		"WM (1)"="#D1C4B0", "WM (2)"="#ede8e0", "GABA"="#9377AC", "low UMI"="grey"),
        precast_k7_HM = c("Vasc"="#911223", "L1"="#cfa45c", "L2"="#5D9940", 
		"L3/4"="#5095CD", "L5"="#ddc94e", "L6"="#E45C5F", "WM"="#D1C4B0"),
        precast_k9_HM = c("Vasc"="#911223", "L1"="#cfa45c", "L2"="#5D9940", 
		"L3/4"="#5095CD", "L5"="#ddc94e", "L6"="#E45C5F", 
		"WM"="#D1C4B0", "GABA"="#9377AC", "L2/3?"="grey"),
	precast_k9_1626 = c("Vasc"="#911223", "L1"="#cfa45c", "L2"="#5D9940",
                "L3/4 (1)"="#5095CD", "L3/4 (2)"="#a3c7e4", "L5"="#ddc94e", "L6"="#E45C5F",
                "WM"="#D1C4B0", "low UMI"="grey"),
	precast_k7_1626 = c("L1"="#cfa45c", "L2"="#5D9940", "L3/4"="#5095CD", 
		"L5"="#ddc94e", "L6"="#E45C5F", "WM"="#D1C4B0", "low UMI"="lightgrey")
)
