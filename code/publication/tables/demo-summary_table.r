library(dplyr)
set.seed(123)

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
  summarise(n=n(),
            age_avg= round(mean(age),1), age_sd= round(sd(age),1), #age_p1sd = age_avg+age_sd, age_m1sd= age_avg-age_sd,
	bmi_avg= round(mean(BMI),1), bmi_sd= round(sd(BMI),1), #bmi_p1sd = bmi_avg+bmi_sd, bmi_m1sd= bmi_avg-bmi_sd,
	rin_avg= round(mean(RIN),2), rin_sd= round(sd(RIN),2), #rin_p1sd = rin_avg+rin_sd, rin_m1sd= rin_avg-rin_sd,
	n_nic= sum(Smoking=="Yes"), n_sud= sum(suicide),
            )

write.csv(t1, "processed-data/publication/Table1_demo-summary_counts-means.csv", row.names=F)

c1 = c("age","BMI","RIN")
c2 = c("Smoking","suicide")

t2 = do.call(rbind, lapply(c1, function(x) {
	demo$plot.me= scale(demo[,x])
	m1 = anova(lm(plot.me ~ condition*sex, data=demo))
	data.frame("var"=x, df.total=m1$Df[[4]], df.dx=m1$Df[[1]], test.dx=signif(m1[["F value"]][[1]],3), p.dx=signif(m1[["Pr(>F)"]][[1]],3),
		df.sex=m1$Df[[2]], test.sex=signif(m1[["F value"]][[2]],3), p.sex=signif(m1[["Pr(>F)"]][[2]],3),
		df.dxsex=m1$Df[[3]], test.dxsex=signif(m1[["F value"]][[3]],3), p.dxsex=signif(m1[["Pr(>F)"]][[3]],3))
}))

demo$group = paste(demo$condition, demo$sex)

t3 = do.call(rbind, lapply(c2, function(x) {
	xsq1 = chisq.test(table(demo[,c("condition",x)]))
	xsq2 = chisq.test(table(demo[,c("sex",x)]))
	xsq3 = chisq.test(table(demo[,c("group",x)]))
	data.frame(var=x, df.total=NA, df.dx=xsq1$parameter, test.dx=signif(xsq1$statistic, 3), p.dx=signif(xsq1$p.value, 3),
		df.sex=xsq2$parameter, test.sex=signif(xsq2$statistic, 3), p.sex=signif(xsq2$p.value, 3),
		df.dxsex=xsq3$parameter, test.dxsex=signif(xsq3$statistic, 3), p.dxsex=signif(xsq3$p.value, 3))
})) 

write.csv(rbind(t2,t3), "processed-data/publication/Table1_demo-summary_stats.csv", row.names=F)

