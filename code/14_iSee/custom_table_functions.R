WT_SUMMARY <- function(se, ri, ci) {
  if (is.null(ri)) {
    ri <- rownames(se)
  } else {
    ri <- unique(unlist(ri))
  }
  if (is.null(ci)) {
    ci <- colnames(se)
  } else {
    ci <- unique(unlist(ci))
  }
  
  df1 = do.call(rbind, lapply(c("domain-SP","domain-CT"), function(annot_name) {
    rdata = as.data.frame(metadata(se)[[paste0(annot_name, "_DE")]][ri,])
    rdata = rdata[,c(1:3, grep("whole", colnames(rdata)))]
    colnames(rdata)[6:ncol(rdata)] = substr(colnames(rdata)[6:ncol(rdata)], start=26, stop=50)
    rdata[,6:ncol(rdata)] = apply(rdata[,6:ncol(rdata)], MARGIN=2, function(x) ifelse(x<.05, "t-test adj. p<.05", "t-test NS"))
    
    colnames(rdata)[4:5] = c("dx-sex_F_stat","dx-sex_F_adjp")
    rdata[,4] = signif(rdata[,4], digits=3)
    rdata[,5] = ifelse(rdata[,5]<.001, format(signif(rdata[,5], digits=3), scientific=T), signif(rdata[,5], digits=3))
    
    rdata$annotation = annot_name
    
    
    rdata[,c("gene_name","annotation","dx-sex_F_stat","dx-sex_F_adjp","F_NTC.MDD","F_NTC.BPD","F_MDD.BPD","M_NTC.MDD","M_NTC.BPD","M_MDD.BPD")]
    
  })) %>% as.data.frame()
  
  rownames(df1) <- NULL
  df1$gene_name = factor(df1$gene_name, levels=ri)
  return(df1[order(df1$gene_name),])
  
}


DR_SUMMARY <- function(se, ri, ci, annot_name = c("domain-SP","domain-CT")) {
  if (is.null(ri)) {
    ri <- rownames(se)
  } else {
    ri <- unique(unlist(ri))
  }
  if (is.null(ci)) {
    ci <- colnames(se)
  } else {
    ci <- unique(unlist(ci))
  }
  
  annot_name = match.arg(annot_name)
  
  rdata = as.data.frame(metadata(se)[[paste0(annot_name, "_DE")]][ri,])
  rdata = rdata[,c(1:3, grep("domain", colnames(rdata)))]
  
  colnames(rdata)[4:5] = c("dx-sex-domain_F_stat","dx-sex-domain_F_adjp")
  rdata[,4] = signif(rdata[,4], digits=3)
  rdata[,5] = ifelse(rdata[,5]<.001, format(signif(rdata[,5], digits=3), scientific=T), signif(rdata[,5], digits=3))
  
  
  domains = unique(sapply(strsplit(colnames(rdata)[6:ncol(rdata)], "_"), function(x) unlist(x)[[4]]))
  df1 <- do.call(rbind, lapply(domains, function(x) {
    r_sub = rdata[,c(1:5, grep(x, colnames(rdata)))]
    colnames(r_sub)[6:ncol(r_sub)] = sapply(strsplit(colnames(r_sub)[6:ncol(r_sub)], "_"), function(x) paste(unlist(x)[5:6], collapse="_"))
    r_sub[,6:ncol(r_sub)] = apply(r_sub[,6:ncol(r_sub)], MARGIN=2, function(x) ifelse(x<.05, "t-test adj. p<.05", "t-test NS"))
    r_sub$domain = x
    return(r_sub)
  }))
  df1$annotation = annot_name
  df1 = df1[,c("gene_name","annotation","dx-sex-domain_F_stat","dx-sex-domain_F_adjp","domain","F_NTC.MDD","F_NTC.BPD","F_MDD.BPD","M_NTC.MDD","M_NTC.BPD","M_MDD.BPD")]
  
  rownames(df1) <- NULL
  df1$gene_name = factor(df1$gene_name, levels=ri)
  return(df1[order(df1$gene_name),])
}

addDEGsToRowData <- function(spe_pseudo) {
  l1 <- lapply(c("domain-SP","domain-CT"), function(annot_name) {
    mdata = as.data.frame(metadata(spe_pseudo)[[paste0(annot_name, "_DE")]])
    list("wt"=mdata[,5]<.05 & rowSums(mdata[, grep("whole.tissue_t", colnames(mdata))]<.05)>0,
         "dr"=mdata[,13]<.05 & rowSums(mdata[,grep("domain.restricted_t", colnames(mdata))]<.05)>0)
  })
  stopifnot(identical(names(l1[[1]]$wt), names(l1[[2]]$wt)))
  wt.deg = factor(paste(l1[[1]]$wt, l1[[2]]$wt), levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
                  labels=c("none","domain-CT","domain-SP","domain-SP & domain-CT"))
  stopifnot(identical(names(l1[[1]]$dr), names(l1[[2]]$dr)))
  dr.deg = factor(paste(l1[[1]]$dr, l1[[2]]$dr), levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
                  labels=c("none","domain-CT","domain-SP","domain-SP & domain-CT"))
  
  rdata = rowData(spe_pseudo)
  stopifnot(identical(names(l1[[1]]$wt), rdata$gene_name))
  rdata$whole.tissue_DEG = wt.deg
  identical(names(l1[[1]]$dr), rdata$gene_name)
  rdata$domain.restricted_DEG = dr.deg
  
  return(rdata)
}
