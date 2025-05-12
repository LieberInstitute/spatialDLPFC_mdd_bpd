getMBvSampleInfo <- function(REDCapFile, demoFile) {
	#REDCap checks and load
	if(length(grep(".csv", REDCapFile))==0) REDCapFile = paste0(REDCapFile,".csv")
	if(length(grep("raw-data", REDCapFile))==0) REDCapFile = paste0("raw-data/sample_info/",REDCapFile)
	if(file.exists(REDCapFile)==FALSE) {
		stop(paste("No file of name",REDCapFile,"found. Check file name or move file into raw-data/sample_info/ dir."))
	} else {REDCap <- read.csv(REDCapFile)}
	
	#demographic checks and load
	if(length(grep(".csv", demoFile))==0) demoFile = paste0(demoFile,".csv")
	if(length(grep("raw-data", demoFile))==0) demoFile = paste0("raw-data/sample_info/",demoFile)
	if(file.exists(demoFile)==FALSE) {
		stop(paste("No file of name",demoFile,"found. Check file name or move file into raw-data/sample_info/ dir."))
	} else {demo <- read.csv(demoFile)}

	#format REDCap
	A1 <- subset(REDCap, select = c("slide", "sample_a1", "project_a1","sample_number1_a1"))
	B1 <- subset(REDCap, select = c("slide", "sample_b1", "project_b1","sample_number1_b1"))
	C1 <- subset(REDCap, select = c("slide", "sample_c1", "project_c1","sample_number1_c1"))
	D1 <- subset(REDCap, select = c("slide", "sample_d1", "project_d1","sample_number1_d1"))
	colnames(A1) <- colnames(B1) <- colnames(C1) <- colnames(D1) <- c("slide", "sample", "project","MBv_sample")
	A1$array <- "A1"
	B1$array <- "B1"
	C1$array <- "C1"
	D1$array <- "D1"
	REDCap_table <- rbind(A1, B1, C1, D1)
	REDCap_table <- REDCap_table[order(REDCap_table$slide), ]

	#filter to mbv
	REDCap_MBv <- REDCap_table[which(REDCap_table$project == "spatialDLPFC_MBv_4100"), ]

	#fix slide 301 array info (see note in raw-data/sample_info/README
	REDCap_MBv[REDCap_MBv[["slide"]]=="V13B23-301","array"] = c("D1","C1","B1","A1")

	#fix "Mbv_034" to "MBv_034"
	REDCap_MBv[grep("b", REDCap_MBv[["MBv_sample"]]),"MBv_sample"] = "MBv_034"

	#fix MBv samples with _DO-NO_SEQ
	v1 = unique(REDCap_MBv[["MBv_sample"]])
	names(v1) = v1
	v2 = sapply(strsplit(v1,"_"), length)
	#v2[v2>2]
	#substr(names(v2[v2>2]), start=0, stop=7)
	REDCap_MBv[["MBv_sample"]] = substr(REDCap_MBv[["MBv_sample"]], start=0, stop=7)
	
	#add sample_id
	REDCap_MBv[["sample_id"]] = paste(REDCap_MBv[["slide"]], REDCap_MBv[["array"]], sep="_")

	#add sequencing tech var -- coverted to code based on methods written by Svitlana
	#Samples MBv_001-008 and MBv_013-016 were sequenced on an S4 NovaSeq 6000 (Illumina) at the SC-TC. 
	#MBv_017-024 were run on a S4 NovaSeq 6000 at the SKCCC. 
	#The remaining samples (MBv_009-012 and MBv_025-120) were sequencing at Psomagen over 9 lanes of a 25B Novaseq X. 
	#MBv_121-128, and were sequenced on one lane of a 25B Novaseq X at Psomagen in a fourth sequencing batch
	REDCap_MBv[["seq"]] = ifelse(REDCap_MBv[["MBv_sample"]] %in% c(paste0("MBv_00",1:8), paste0("MBv_0",13:16)), "SC-TC", "other")
	REDCap_MBv[["seq"]] = ifelse(REDCap_MBv[["MBv_sample"]] %in% paste0("MBv_0",17:24), "SKCCC", REDCap_MBv[["seq"]])
	REDCap_MBv[["seq"]] = ifelse(REDCap_MBv[["MBv_sample"]] %in% c("MBv_009", paste0("MBv_0",10:12), paste0("MBv_0",25:99), paste0("MBv_",100:120)), 
				"Psomagen-1", REDCap_MBv[["seq"]] )
	REDCap_MBv[["seq"]] = ifelse(REDCap_MBv[["MBv_sample"]] %in% paste0("MBv_",121:128), "Psomagen-2", REDCap_MBv[["seq"]])
	
	#add some extra columns
	r1 = c("V13F27-338", "V13F27-348", "V13Y10-020", "V13Y10-021", "V13Y10-022", "V13Y10-023")
	r3 = c("V13B23-282","V13B23-279") #V13Y10-020 --> V13B23-282, V13B23-331 --> V13B23-279
	#r2_1 = "V13B23-283" #this slide was actually sequenced and sent out in the same batch
	#r2 = setdiff(outDF$slide, c(r1, r3, r2_1))
	r2 = setdiff(REDCap_MBv[["slide"]], c(r1, r3))
	REDCap_MBv[["round"]] = factor(REDCap_MBv[["slide"]], levels=c(r1, r2, r3),
		#levels=c(r1, r2, r2_1, r3), 
		labels=c(rep("r1", length(r1)),
			rep("r2", length(r2)),
			#"r2_1",
			rep("r3", length(r3))
		)
	)

	#combine with demo data
        outDF = merge(demo, REDCap_MBv[,c("sample","sample_id","slide","array","MBv_sample","seq","round")],
		by.x="brnum", by.y="sample")

	return(outDF)
}
