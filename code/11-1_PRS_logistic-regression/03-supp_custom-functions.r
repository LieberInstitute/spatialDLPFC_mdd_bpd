prs_wrapper <- function(risk_set) {
  # load PRS
  prs.df = read.csv(paste0("raw-data/PRS/PRS_",risk_set,".csv"))
  colnames(prs.df)[1] = "brnum"
  prs.df$risk = risk_set
  
  prs.df = prs.df[,c("brnum","risk",paste0("p.cutoff.",c("1e.08","1e.07","1e.06","1e.05")))]
  
  cutoffs = grep("p\\.cutoff", colnames(prs.df), value=T)
  
  outList = lapply(cutoffs, function(x) {
    colnames(prs.df)[colnames(prs.df)==x] = "PRS"
    prs.df$PRS_scaled = scale(prs.df$PRS)
    return(prs_processing(prs.df))
  })
  names(outList) <- cutoffs
  return(outList)
}

prs_processing <- function(prs.df) {
  # load in donor metadata
  cdata = read.csv("raw-data/sample_info/DLPFC_cross-disorders_demographics_MBv.csv")
  cdata$RIN_scaled = scale(cdata$RIN)
  
  prs.df = left_join(prs.df, cdata) %>%
    mutate(sex=factor(sex, levels=c("F","M")),
           condition=factor(condition, levels=c("NTC","MDD","BPD")))
  
  # pull SNP PCs
  t1 = read.table("/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/11_eQTL_coloc/seurat/tqtl_in/astro.gene.covars.txt", row.names=1)
  
  t2 = as.data.frame(apply(t(t1[-1,]), MARGIN=2, as.numeric))
  t2$brnum = as.character(t1[1,])
  
  snp.df = t2[,c("brnum",grep("snp", colnames(t2), value=T))]
  snp.df2 = apply(snp.df[,-1], MARGIN=2, scale)
  colnames(snp.df2) = paste0(colnames(snp.df2),"_scaled")
  snp.df = cbind(snp.df, snp.df2)
  stopifnot(nrow(snp.df)==119)
  
  # merge snp covars
  all.df = left_join(prs.df, snp.df)
  
  risk_set = as.character(factor(unique(prs.df$risk), levels=c("MDD","Bipolar","SCZ"), labels=c("MDD","BPD","SCZ")))
  
  ol = list()
  
  set.seed(123)
  cond.df = filter(all.df, condition %in% c("NTC","MDD")) %>% mutate(condition= droplevels(condition))
  ol[[length(ol)+1]] = glm(condition ~ PRS_scaled + age + sex + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled, 
                           data=cond.df, family=binomial)
  names(ol)[length(ol)] = "dxMDD"
    
  set.seed(123)
  cond.df = filter(all.df, condition %in% c("NTC","BPD")) %>% mutate(condition= droplevels(condition))
  ol[[length(ol)+1]] = glm(condition ~ PRS_scaled + age + sex  + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled, 
                           data=cond.df, family=binomial)
  names(ol)[length(ol)] = "dxBPD"
  
  set.seed(123)
  all.df2 = mutate(all.df, condition= factor(condition, levels=c("NTC","MDD","BPD"), labels=c("NTC","DX","DX")))
  ol[[length(ol)+1]] = glm(condition ~ PRS_scaled + age + sex + snpPC1_scaled + snpPC2_scaled + snpPC3_scaled + snpPC4_scaled + snpPC5_scaled, 
                           data=all.df2, family=binomial)
  names(ol)[length(ol)] = "dxAny"
  
  return(ol)
}


sigCoef <- function(modelList) {
  lapply(modelList, function(y) {
    z = coef(summary(y))
    is_prs_sig = z["PRS_scaled", "Pr(>|z|)"]<.05
    if(is_prs_sig) {
      return(exp(coef(y))[["PRS_scaled"]])
    } else return("NS")
  })
}

getPval <- function(modelList) {
  unlist(lapply(modelList, function(y) {
    z = coef(summary(y))
    return( z["PRS_scaled", "Pr(>|z|)"])
  }))
}

meltPvalDF <- function(modelList) {
  pmtx = sapply(modelList, getPval)
  padjmtx = apply(pmtx, MARGIN=1, FUN=p.adjust, method="BH")
  df1 = tibble::rownames_to_column(as.data.frame(t(pmtx)), var="cutoff") %>%
    tidyr::pivot_longer(c("dxMDD","dxBPD","dxAny"), names_to="dx", values_to="raw_p")
  df2 = tibble::rownames_to_column(as.data.frame(padjmtx), var="cutoff") %>%
    tidyr::pivot_longer(c("dxMDD","dxBPD","dxAny"), names_to="dx", values_to="adj_p")
  merge(df1, df2)
}

getConfInt <- function(modelList) {
  df1 = do.call(rbind, lapply(modelList, function(model1) {
    suppressMessages({
      data.frame("OR"=exp(coef(model1))[["PRS_scaled"]],
                 "CI05"=exp(confint(model1))["PRS_scaled",1],
                 "CI95"=exp(confint(model1))["PRS_scaled",2])
    })
  }))
  df1$dx = rownames(df1)
  rownames(df1) <- NULL
  return(df1[,c(4,1:3)])
}
