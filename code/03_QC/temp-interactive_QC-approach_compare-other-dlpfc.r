library(SpatialExperiment)
library(dplyr)
library(ggplot2)
#my raw data number of in tissue low umi spots
cdata = read.csv("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
dim(cdata)
#599034     38
table(cdata$in_tissue)
#FALSE   TRUE 
#46072 552962
table(cdata[cdata$in_tissue,"lowumi"])
#FALSE   TRUE
#542562  10400
#10400/(10400+542562) #0.0188

filter(cdata, remove_spots==FALSE) %>% group_by(sample_id) %>% 
  summarise(med.umi=median(sum_umi)) %>% pull(med.umi) %>% sort()
### 5 out of 120 sections median sum_umi <1000 # 0.0417
### 15 out of 120 sections median sum_umi <1500 # 0.125 or 1 out of 8
#[1]  577.0  676.0  850.0  938.0  957.5 1151.0 1244.0 1271.0 1289.0 1291.0
#[11] 1309.0 1357.5 1367.0 1443.0 1462.0 1512.0 1671.0 1705.5 1738.5 1760.5
#[21] 1796.0 1848.5 1898.5 1915.0 1924.0 1935.0 1936.0 1962.0 1973.0 1991.0

filter(cdata, remove_spots==FALSE) %>% group_by(sample_id) %>% 
  summarise(med.genes=median(sum_gene)) %>% pull(med.genes) %>% sort()
### 15 out of 120 sections median sum_gene <1000 # 0.125 or 1 out of 8
### 19 out of 120 sections median sum_gene >2000 # 0.158
#[1]  375.0  476.0  561.0  629.0  669.5  698.0  721.0  751.0  799.0  806.0
#[11]  820.0  856.0  866.5  937.0  943.0 1001.0 1002.0 1021.0 1025.0 1048.5
#[21] 1095.0 1113.0 1114.0 1124.0 1170.0 1171.5 1193.0 1196.0 1223.0 1238.0
#...
#[101] 1950.0 2040.0 2050.0 2057.0 2081.5 2109.5 2132.0 2139.5 2150.0 2160.0
#[111] 2201.5 2210.0 2250.0 2357.0 2416.0 2439.5 2444.0 2456.5 2483.0 3272.0


#########################################################################################
#########################################################################################

"/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100" #Boyi's dataset
#https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/02_visium_qc/assemble_qc_spe.r
spe_raw <- readRDS("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/01_build_spe/raw_spe_wo_SPG_N63.rds")
#having memory issues with 30GB session so keeping just colData
cdata_raw = as.data.frame(colData(spe_raw))
rm(spe_raw)

dim(spe_raw)
#36601 314495
table(spe_raw$in_tissue)
#FALSE   TRUE 
#31688 282807 

#number of in tissue low umi spots
table(colData(spe_raw)[spe_raw$in_tissue,"sum_umi"]<=100)
#FALSE   TRUE 
#282256    551 
# 551/(551+282256) #0.00195
table(colData(spe_raw)[spe_raw$in_tissue,"sum_umi"]<=200)
#FALSE   TRUE 
#282018    789
# 789/(789+282018) #0.00279


filter(as.data.frame(colData(spe_raw)), in_tissue==TRUE) %>%
  group_by(sample_id) %>% summarise(med.umi=median(sum_umi)) %>% pull(med.umi) %>% sort()
### even raw (in_tissue filters) 0 out of 63 sections with median umi <1500
### even raw 43 out of 63 sections with median umi >2000
#[1] 1579.5 1710.5 1786.0 1913.5 2153.0 2191.0 2328.0 2382.0 2474.0 2523.0

#https://github.com/LieberInstitute/spatialDLPFC_SCZ/blob/main/code/analysis/02_visium_qc/post_QC_stat.r
spe_pnn <- readRDS("/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/rds/02_visium_qc/qc_spe_wo_spg_N63.rds")
#having memory issues with 30GB session so keeping just colData
cdata_pnn = as.data.frame(colData(spe_pnn))
rm(spe_pnn)

dim(spe_pnn)
#32540 279806
table(spe_pnn$in_tissue)
#TRUE 
#279806 

#number of in tissue low umi spots
table(colData(spe_pnn)$sum_umi<=100)
table(colData(spe_pnn)$sum_umi<=200)


filter(as.data.frame(colData(spe_raw)), key %in% spe_pnn$key) %>% dim() #279806
filter(as.data.frame(colData(spe_raw)), key %in% spe_pnn$key) %>%
  group_by(sample_id) %>% summarise(med.umi=median(sum_umi)) %>% pull(med.umi) %>% sort()
### 0 out of 63 samples with less than median umi <1500
#[1] 1589.0 1722.5 1793.0 1916.0 2154.5 2233.0 2329.0 2422.0 2516.0 2523.0

filter(as.data.frame(colData(spe_raw)), key %in% spe_pnn$key) %>%
  group_by(sample_id) %>% summarise(med.genes=median(sum_gene)) %>% pull(med.genes) %>% sort()
### 2 out of 63 sections with median sum_gene <1000
### 20 out of 63 sections with median sum_gene >2000
#[1]  907.0  981.0 1074.0 1180.0 1199.0 1250.0 1353.0 1356.0 1384.0 1403.0
#...
#[41] 1957.0 1991.5 2010.0 2057.0 2083.5 2093.0 2100.0 2129.0 2163.0 2214.0
#[51] 2268.0 2353.0 2396.0 2419.0 2503.5 2596.0 2602.0 2616.0 2934.0 2961.0
#[61] 2983.5 3321.0 3394.0


intersect(cdata$brnum, spe_pnn$brnum)
#[1] "Br1092" "Br1113" "Br1204" "Br2720" "Br5228" "Br5276" "Br5314" "Br5367"
#[9] "Br5395" "Br5400" "Br5436" "Br5472" "Br5599" "Br5622" "Br5639" "Br6297"
#[17] "Br6432" "Br8042" "Br8218" "Br8325" "Br8433" "Br8492" "Br8667"

pnn.summary = right_join(as.data.frame(colData(spe_raw)), as.data.frame(colData(spe_pnn)), by="key") %>%
  filter(brnum %in% intersect(cdata$brnum, spe_pnn$brnum)) %>%
  group_by(brnum) %>% summarise(med.umi_pnn=median(sum_umi))
mbv.summary = filter(cdata, brnum %in% intersect(cdata$brnum, spe_pnn$brnum), in_tissue==TRUE) %>% 
  group_by(brnum) %>% mutate(med.umi_unfilt=median(sum_umi)) %>%
  filter(remove_spots==FALSE) %>% group_by(brnum, sample_id, med.umi_unfilt) %>%
  summarise(med.umi_filt=median(sum_umi))

comb.summary = left_join(pnn.summary, mbv.summary[,c("brnum","med.umi_unfilt","med.umi_filt","sample_id")], by="brnum")
print(comb.summary, n=23)
#brnum  med.umi_pnn med.umi_unfilt med.umi_filt sample_id    
#1 Br1092       4305           3677         3684  V13B23-380_C1
#2 Br1113       4431           2451         2454. V13B23-342_A1
#3 Br1204       4090.          2530         2589  V13B23-302_C1
#4 Br2720       2759           3005         3010  V13B23-381_D1
#5 Br5228       3826.          3150         3152  V13B23-302_D1
#6 Br5276       3953           2938.        2941  V13B23-332_B1
#7 Br5314       3068           2656         2812  V13F27-338_B1
#8 Br5367       7446           4799         4811  V13B23-380_D1
#9 Br5395       4056.          3903         3912. V13B23-311_A1
#10 Br5400       5101           2356.        2364. V13Y10-022_D1
#11 Br5436       5634.          5017         5217  V13B23-334_D1
#12 Br5472       6447           3052         3140. V13B23-310_A1
#13 Br5599       2329           1872         1936  V13Y10-022_C1
#14 Br5622       3014           2457         2461  V13B23-310_B1
#15 Br5639       3811           1216         1244  V13Y10-021_C1
#16 Br6297       3792.          2673         2686  V13B23-311_B1
#17 Br6432       3076           5339         5340  V13B23-308_D1
#18 Br8042       4139           2884         3224  V13B23-330_B1
#19 Br8218       2516           1268.        1271  V13B23-340_B1
#20 Br8325       3775           1755         1760. V13B23-343_C1
#21 Br8433       2154.          2187         2221  V13B23-279_B1
#22 Br8492       2597           2488         2488. V13B23-327_C1
#23 Br8667       2685           4322         4326  V13B23-332_A1

#########################################################################################
#########################################################################################

"/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC" #Louise's dataset
#https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/01_build_spe/01_build_spe.R
#save(spe_raw, file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_raw_final.Rdata"))
#save(spe, file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_final.Rdata"))
#save(spe, file = here::here("processed-data", "rdata", "spe", "01_build_spe", "spe_filtered_final.Rdata"))
library(SpatialExperiment)
library(dplyr)

### RAW DATA
load("processed-data/rdata/spe/01_build_spe/spe_raw_final.Rdata")
dim(spe_raw)
#36601 149757
table(spe_raw$in_tissue)
#FALSE   TRUE 
#30957 118800 

#number of in tissue low umi spots
table(colData(spe_raw)[spe_raw$in_tissue,"sum_umi"]<=100)
#FALSE   TRUE 
#117054   1746
# 1746/(1746+117054) #0.0147

### FINAL 1 (looks like just removed off tissue spots removed plus 7 more spots)
load("processed-data/rdata/spe/01_build_spe/spe_final.Rdata")
dim(spe)
#28916 118793
table(spe$in_tissue)
# TRUE 
#118793 

table(spe$sum_umi<=100)
#FALSE   TRUE 
#117054   1739
#1739/(1739+117054) #0.0146

### FILTERED QC
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final.Rdata")
#having memory issues with 30GB session so keeping just colData
cdata_dlpfc = as.data.frame(colData(spe))
rm(spe)

dim(spe)
#28916 113927
table(spe$in_tissue)
# TRUE 
#113927

table(spe$sum_umi<=100)
#FALSE   TRUE 
#113739    188
# 188/(188+113739) #0.00165

table(spe$sum_umi<=200)
#FALSE   TRUE 
#112782   1145 
# 1145/(1145+112782) #0.0101

#my sum_umi<=200 after my proposed QC filters
table(cdata[cdata$remove_spots==FALSE,"sum_umi"]<=200)
#FALSE   TRUE 
#532285   7123
# 7123/(7123+532285) #0.0132


group_by(as.data.frame(colData(spe)), sample_id) %>%
  summarise(med.umi=median(sum_umi)) %>% pull(med.umi) %>% sort()
### 1 out of 30 sections median sum_umi <1000  # 0.033
### 1 out of 30 sections median sum_umi <1500  # 0.033
#[1]  347.0 1626.0 1648.5 1649.0 1741.0 1788.0 1925.0 1933.0 2035.5 2080.5
#[11] 2147.0 2184.0 2197.0 2233.0 2328.0 2351.0 2438.0 2462.0 2535.0 2599.0
#[21] 2812.0 2848.5 3252.0 3468.0 3498.5 3874.0 4152.5 4179.0 4284.0 4765.0


### CHECK COMPARISON OF SAME BRAIN DONORS
intersect(cdata$brnum, spe$subject)
#[1] "Br2720" "Br6432" "Br6471" "Br8325" "Br8492" "Br8667"

filter(as.data.frame(colData(spe)), subject %in% intersect(cdata$brnum, spe$subject)) %>% 
  mutate(lowumi=sum_umi<=100) %>%
  group_by(subject) %>% mutate(med.umi_subject=median(sum_umi)) %>%
  group_by(subject, sample_id, med.umi_subject) %>%
  summarise(med.umi_sample=median(sum_umi), n_lowumi=sum(lowumi))

#subject sample_id   med.umi_subject med.umi_sample n_lowumi
#1 Br2720  Br2720_ant             1543           347       188
#2 Br2720  Br2720_mid             1543          2036.        0
#3 Br2720  Br2720_post            1543          2535         0
#4 Br6432  Br6432_ant             2957          1741         0
#5 Br6432  Br6432_mid             2957          3874         0
#6 Br6432  Br6432_post            2957          4765         0
#7 Br6471  Br6471_ant             2056          1788         0
#8 Br6471  Br6471_mid             2056          2462         0
#9 Br6471  Br6471_post            2056          1933         0
#10 Br8325  Br8325_ant             1927          2184         0
#11 Br8325  Br8325_mid             1927          2080.        0
#12 Br8325  Br8325_post            1927          1648.        0
#13 Br8492  Br8492_ant             1909          1925         0
#14 Br8492  Br8492_mid             1909          2233         0
#15 Br8492  Br8492_post            1909          1626         0
#16 Br8667  Br8667_ant             2939          2438         0
#17 Br8667  Br8667_mid             2939          4179         0
#18 Br8667  Br8667_post            2939          2599         0

filter(cdata, brnum %in% intersect(cdata$brnum, spe$subject), in_tissue==TRUE) %>% 
  group_by(brnum) %>% mutate(med.umi_unfilt=median(sum_umi)) %>%
  filter(remove_spots==FALSE) %>% group_by(brnum, sample_id, med.umi_unfilt) %>%
  summarise(med.umi_subject=median(sum_umi))
#brnum  sample_id     med.umi_unfilt med.umi_subject
#1 Br2720 V13B23-381_D1          3005            3010 #closest to post section (but higher still)
#2 Br6432 V13B23-308_D1          5339            5340 #closest to post section (but higher still)
#3 Br6471 V13B23-381_C1          2070.           2078 #closest to post section
#4 Br8325 V13B23-343_C1          1755            1760 #closest to post section
#5 Br8492 V13B23-327_C1          2488            2488 #closest to mid section (but higher still)
#6 Br8667 V13B23-332_A1          4322            4326 #closest to mid section (but higher still)



group_by(as.data.frame(colData(spe)), sample_id) %>%
  summarise(med.genes=median(sum_gene)) %>% pull(med.genes) %>% sort()
### 1 out of 30 sections median sum_gene <1000 # 0.033
### 7 out of 30 sections median sum_gene >2000 # 0.233
#[1]  262.0 1062.5 1065.0 1114.0 1142.0 1151.0 1207.5 1233.0 1278.0 1280.0
#[11] 1320.0 1354.0 1367.0 1393.0 1408.0 1446.0 1448.5 1485.0 1510.0 1564.0
#[21] 1662.0 1692.5 1719.0 2048.0 2101.0 2160.0 2164.0 2177.5 2451.0 2546.0

#########################################################################################
#########################################################################################

### COMPARE PROJECTS
comb.df = bind_rows(group_by(as.data.frame(colData(spe)), sample_id, region) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      rename(group=region) %>% mutate(source="spatialDLPFC"),
                    filter(cdata, in_tissue==TRUE) %>% group_by(sample_id) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="MBv", group="unfiltered"),
                    filter(cdata, remove_spots==FALSE) %>% group_by(sample_id) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="MBv", group="filtered"),
                    filter(as.data.frame(colData(spe_raw)), in_tissue==TRUE) %>% group_by(sample_id) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="PNN/SCZ", group="unfiltered"),
                    filter(as.data.frame(colData(spe_raw)), key %in% spe_pnn$key) %>% group_by(sample_id) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="PNN/SCZ", group="filtered")
) %>% mutate(group=factor(group, levels=c("unfiltered","filtered","anterior","middle","posterior")))

ggplot(comb.df, aes(x=group, y=med.umi))+
  geom_point(data=filter(comb.df, source=="spatialDLPFC"))+
  ggbeeswarm::geom_quasirandom(data=filter(comb.df, source!="spatialDLPFC"))+
  facet_wrap(vars(source), scales="free_x")+
  labs(y="sample median sum_umi", title="Library size")+
  theme(text=element_text(size=16))

ggplot(comb.df, aes(x=group, y=med.genes))+
  geom_point(data=filter(comb.df, source=="spatialDLPFC"))+
  ggbeeswarm::geom_quasirandom(data=filter(comb.df, source!="spatialDLPFC"))+
  facet_wrap(vars(source), scales="free_x")+ylim(0,4000)+
  labs(y="sample median sum_gene", title="Detected genes")+
  theme(text=element_text(size=16))

### DONOR LEVEL MATCH
all.donors = intersect(intersect(cdata$brnum, cdata_dlpfc$subject), cdata_pnn$brnum)
all.comb.df = bind_rows(filter(cdata_dlpfc, subject %in% all.donors) %>% group_by(subject, region) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      rename(group=region, brnum=subject) %>% mutate(source="spatialDLPFC") %>% ungroup(),
                    filter(cdata, in_tissue==TRUE, brnum %in% all.donors) %>% group_by(brnum) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="MBv", group="unfiltered")  %>% ungroup(),
                    filter(cdata, remove_spots==FALSE,  brnum %in% all.donors) %>% group_by(brnum) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="MBv", group="filtered")  %>% ungroup(),
                    left_join(cdata_raw, distinct(cdata_pnn[,c("sample_id","brnum")]), by="sample_id") %>% 
                      filter(in_tissue==TRUE, brnum %in% all.donors) %>% group_by(brnum) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="PNN/SCZ", group="unfiltered")  %>% ungroup(),
                    right_join(cdata_raw, cdata_pnn[,c("key","brnum")], by="key") %>% filter(brnum %in% all.donors) %>%
                      group_by(brnum) %>% summarise(med.umi=median(sum_umi), med.genes=median(sum_gene)) %>%
                      mutate(source="PNN/SCZ", group="filtered") %>% ungroup()
) %>% mutate(group=factor(group, levels=c("unfiltered","filtered","anterior","middle","posterior")))

ggplot(all.comb.df,  aes(x=group, y=med.umi, color=source))+
  geom_point(size=3)+facet_wrap(vars(brnum))+ylim(0,6000)+
  labs(y="median sum_umi", title="Library size (donor comparison)")+
  theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

ggplot(all.comb.df,  aes(x=group, y=med.genes, color=source))+
  geom_point(size=3)+facet_wrap(vars(brnum))+ylim(0,3000)+
  labs(y="median sum_gene", title="Detected genes (donor comparison)")+
  theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

