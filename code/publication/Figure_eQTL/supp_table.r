library(dplyr)

coloc = readxl::read_excel("processed-data/11_eQTL_coloc/seurat/final/coloc_results.xlsx", sheet = "coloc_pass") %>%
  mutate(disorder= factor(disorder, levels=c("MDD","BD","SCZD")))


col1 = distinct(coloc, disorder, gene_name, DEG, lead_variant_gwas_strict, lead_variant_gwas_exploratory, gene_gwas_list)

filter(col1, lead_variant_gwas_strict==T) %>% group_by(DEG, disorder) %>% tally()
filter(col1, lead_variant_gwas_exploratory==T) %>% group_by(DEG, disorder) %>% tally()
filter(col1, gene_gwas_list==T) %>% group_by(DEG, disorder) %>% tally()

filter(coloc, lead_variant_gwas_strict==T, DEG==1, disorder=="SCZD") %>%
  select(context, gene_name, lead_snp, lead_rsid)
filter(coloc, lead_variant_gwas_exploratory==T, DEG==1, disorder=="MDD") %>%
  select(context, gene_name, lead_snp, lead_rsid)
filter(coloc, lead_variant_gwas_exploratory==T, DEG==1, disorder=="BD") %>%
  select(context, gene_name, lead_rsid)


filter(coloc, lead_variant_gwas_exploratory==T | lead_variant_gwas_strict==T, disorder=="MDD") %>%
  group_by(lead_snp, lead_rsid) %>% summarize(exploratory=sum(lead_variant_gwas_exploratory), strict=sum(lead_variant_gwas_strict))
# nrow = 19 different variants at either threshold

filter(coloc, lead_variant_gwas_exploratory==T | lead_variant_gwas_strict==T, disorder=="BD") %>%
  group_by(lead_snp, lead_rsid) %>% summarize(exploratory=sum(lead_variant_gwas_exploratory), strict=sum(lead_variant_gwas_strict))
# nrow = 14 different variants at either threshold
