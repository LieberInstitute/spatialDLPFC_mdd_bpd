library(dplyr)

demo = read.csv("processed-data/publication/supp_tables/demographics.csv")

df1 = read.csv("raw-data/sample_info/MBv_suicide.csv")
table(df1$PrimaryDx)

# 6 donors not included in study
filter(df1, !BrNum %in% demo$brnum)
filter(df1, !BrNum %in% demo$brnum)$BrNum

demo = left_join(demo, df1[,c("BrNum","Manner.Of.Death","Violent.method.")], by=c("brnum"="BrNum"))

demo$suicide = !is.na(demo$Manner.Of.Death)

table(demo$suicide)
table(demo$Smoking)


t1 = mutate(demo, condition=factor(condition, levels=c("NTC","MDD","BD"))) %>%
  group_by(condition, sex) %>%
  summarise(n=n(), n_nic= sum(Smoking=="Yes"), n_sud= sum(suicide),
            age_avg= round(mean(age),3), age_sd= round(sd(age),3), age_p1sd = age_avg+age_sd, age_m1sd= age_avg-age_sd,
            )

write.csv(t1, "processed-data/publication/Table1_demo-summary.csv", row.names=F)

