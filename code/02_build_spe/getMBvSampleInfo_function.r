getMBvSampleInfo <- function(REDCapFile, demoFile) {
	#REDCap checks and load
	if(length(grep(".csv", REDCapFile))==0) REDCapFile = paste0(REDCapFile,".csv")
	if(length(grep("raw-data", REDCapFile))==0) REDCapFile = paste0("raw-data/sample_info/",REDCapFile)
	if(file.exists(REDCapFile)==FALSE) {
		stop(paste("No file of name",REDCapFile,"found. Check file name or move file into raw-data/sample_info dir."))
	} else {REDCap <- read.csv(REDCapFile)}
	
	#demographic checks and load
	if(length(grep(".csv", demoFile))==0) demoFile = paste0(demoFile,".csv")
	if(length(grep("raw-data", demoFile))==0) demoFile = paste0("raw-data/sample_info/",demoFile)
	if(file.exists(demoFile)==FALSE) {
		stop(paste("No file of name",demoFile,"found. Check file name or move file into raw-data/sample_info dir."))
	} else {demo <- read.csv(demoFile)}

	#format REDCap
	A1 <- subset(REDCap, select = c("slide", "sample_a1", "project_a1"))
	B1 <- subset(REDCap, select = c("slide", "sample_b1", "project_b1"))
	C1 <- subset(REDCap, select = c("slide", "sample_c1", "project_c1"))
	D1 <- subset(REDCap, select = c("slide", "sample_d1", "project_d1"))
	colnames(A1) <- colnames(B1) <- colnames(C1) <- colnames(D1) <- c("slide", "sample", "project")
	A1$array <- "A1"
	B1$array <- "B1"
	C1$array <- "C1"
	D1$array <- "D1"
	REDCap_table <- rbind(A1, B1, C1, D1)
	REDCap_table <- REDCap_table[order(REDCap_table$slide), ]

	#filter to mbv
	REDCap_MBv <- REDCap_table[which(REDCap_table$project == "spatialDLPFC_MBv_4100"), ]

	#combine with demo data
	outDF = merge(demo, REDCap_MBv[,c("sample","slide","array")], by.x="brnum", by.y="sample")
	
	#add some extra columns
	outDF$sample_id = paste(outDF$slide, outDF$array, sep="_")
	r1 = c("V13F27-338", "V13F27-348", "V13Y10-020", "V13Y10-021", "V13Y10-022", "V13Y10-023")
	r3 = c("V13B23-282","V13B23-279") #V13Y10-020 --> V13B23-282, V13B23-331 --> V13B23-279
	r2_1 = "V13B23-283"
	r2 = setdiff(outDF$slide, c(r1, r3, r2_1))
	outDF$round = factor(outDF$slide, levels=c(r1, r2, r2_1, r3), 
		labels=c(rep("r1", length(r1)),
			rep("r2", length(r2)),
			"r2_1",
			rep("r3", length(r3))
		)
	)
	return(outDF)
}
