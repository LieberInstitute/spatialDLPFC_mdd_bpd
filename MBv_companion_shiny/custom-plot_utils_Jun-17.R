suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(rstatix)
  library(ggpubr)
  library(gridExtra)
})

# whole tissue dataframe
extractLogcounts <- function(gene_name, spe_pseudo, annot_name) {
  spe_pseudo = spe_pseudo[,spe_pseudo$annotation==annot_name]
  spe_pseudo[[gene_name]] = logcounts(spe_pseudo)[gene_name,]
  df = as.data.frame(colData(spe_pseudo)[,c("sample_id","condition","sex","domain",gene_name)])
  df$domain = droplevels(df$domain)
  if(annot_name=="domain-SP") df$domain = factor(df$domain, levels=c("L1","L2","L3.4","L5","L6","WM"))
  if(annot_name=="domain-CT") df$domain = factor(df$domain, levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))
  return(df)
}


# summarize for crossbar
crossbarLogcounts <- function(logcount_DF) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  df = group_by(logcount_DF, condition, sex, domain) %>% 
    summarise(yavg=mean(plot.gene), ysd=sd(plot.gene), n=n(), y_min= yavg-ysd, y_max= yavg+ysd) 
  return(df)
}


# create domain-restricted ggplot base
plotDomainRestricted <- function(logcount_DF, summary_DF, annot_name, spe_pseudo, ymax1) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  ## get limits (just ceiling of max for now, will need to change later with ggpubr)
  #ymax1 = ceiling(max(logcount_DF$plot.gene))
  # if mean-SD is less than zero, extend lower limit of y-axis
  if(min(summary_DF$y_min)<0) {
    ymin1 = floor(min(summary_DF$y_min))
  } else {
    ymin1 = 0
  }
  
  
  # subtitle formating based on sig F test
  subtitle_face = ifelse(isFsig(save_name, logcount_DF, spe_pseudo, annot_name),"bold","plain")
  
  dx.pal = c("NTC"="#7f7f7f","MDD"="#FB8861","BPD"="#9260b2")
  
  if(annot_name=="domain-SP") {
    x_labels = c("L1","L2","L3/4","L5","L6","WM")
  } 
  if(annot_name=="domain-CT") {
    x_labels = c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")
  }
  
  ggplot(logcount_DF, aes(x=domain))+
    geom_quasirandom(aes(y=plot.gene, color=condition), dodge.width=.8, cex=1)+
    scale_color_manual("DX", values=dx.pal, labels=c("NTC"="NTC","MDD"="MDD","BPD"="BD"))+
    geom_crossbar(data=summary_DF, aes(y=yavg, ymin=y_min, ymax=y_max, group=condition),
                  color="black", fill="transparent", position=position_dodge(width=.8), width=.6)+
    facet_grid(cols=vars(sex), labeller= as_labeller(c("F"="Female","M"="Male")))+
    ylim(ymin1,ymax1)+scale_x_discrete(labels=x_labels)+
    labs(title=save_name, subtitle=paste(annot_name, getFname(save_name, logcount_DF, spe_pseudo, annot_name), sep=" "), 
         y="logcounts", x=annot_name)+
    theme_bw()+theme(plot.title=element_text(face="italic"), plot.subtitle = element_text(face=subtitle_face),
                          text=element_text(size=16))
}

# whole tissue ggplot base 
plotWholeTissue <- function(logcount_DF, summary_DF, ymax1, ...) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  ## get limits (just ceiling of max for now, will need to change later with ggpubr)
  #ymax1 = ceiling(max(logcount_DF$plot.gene))
  # if mean-SD is less than zero, extend lower limit of y-axis
  if(min(summary_DF$y_min)<0) {
    ymin1 = floor(min(summary_DF$y_min))
  } else {
    ymin1 = 0
  }
  
  dx.pal = c("NTC"="#7f7f7f","MDD"="#FB8861","BPD"="#9260b2")
  
  ggplot(logcount_DF, aes(x=sex))+
    geom_quasirandom(aes(y=plot.gene, color=condition), dodge.width=.8, cex=1)+
    scale_x_discrete(labels=c("Female","Male"))+
    scale_color_manual("DX", values=dx.pal, labels=c("NTC"="NTC","MDD"="MDD","BPD"="BD"))+
    geom_crossbar(data=summary_DF, aes(y=yavg, ymin=y_min, ymax=y_max, group=condition),
                  color="black", fill="transparent", position=position_dodge(width=.8), width=.6)+
    facet_grid(cols=vars(domain), 
               labeller= as_labeller(c("all-SP"=paste0("domain-SP\n", getFname(save_name, filter(logcount_DF, domain=="all-SP"), ..., "domain-SP")),
                                       "all-CT"=paste0("domain-CT\n", getFname(save_name, filter(logcount_DF, domain=="all-CT"), ..., "domain-CT")))))+
    ylim(ymin1,ymax1)+
    labs(title=save_name, 
         y="logcounts", x="")+
    theme_bw()+theme(plot.title=element_text(face="italic"), text=element_text(size=16))
}


# check if gene was run in model
inModel <- function(gene_name, spe_pseudo, annot_name) {
  gene_name %in% rownames(metadata(spe_pseudo)[[paste0(annot_name, "_DE")]])
}

# extract F statistics for labeling
getFname <- function(gene_name, logcount_DF, spe_pseudo, annot_name) {
  if(!inModel(gene_name, spe_pseudo, annot_name)) return("(gene not run in model, excluded by edgeR filters)") 
  rdata= metadata(spe_pseudo)[[paste0(annot_name, "_DE")]][gene_name,]
  
  if(length(unique(logcount_DF$domain))==1) {
    fstat = rdata$whole.tissue_F_stat
    fpadj = rdata$whole.tissue_F_adj.P.Val
  } else {
    fstat = rdata$domain.restricted_F_stat
    fpadj = rdata$domain.restricted_F_adj.P.Val
  }
  
  #format stats for title
  if(fpadj<.001) {
    format_p = format(fpadj, scientific=T, digits=2)
  } else {
    format_p = signif(fpadj, 2)
  }
  
  format_f = signif(fstat, 3)
  
  #return(paste0(gene_name, " (F= ", format_f,", adj. p= ", format_p, ")"))
  return(paste0("(F= ", format_f,", adj. p= ", format_p, ")"))
}

# conditional F sig for formatting
isFsig <- function(gene_name, logcount_DF, spe_pseudo, annot_name) {
  if(!inModel(gene_name, spe_pseudo, annot_name)) return(FALSE)
  rdata= metadata(spe_pseudo)[[paste0(annot_name, "_DE")]][gene_name,]
  
  if(length(unique(logcount_DF$domain))==1) {
    fpadj = rdata$whole.tissue_F_adj.P.Val
  } else {
    fpadj = rdata$domain.restricted_F_adj.P.Val
  }
  
  return(fpadj<.05)
}

# extract t-test results
getTstats <- function(gene_name, spe_pseudo, annot_name) {
  rdata= as.data.frame(metadata(spe_pseudo)[[paste0(annot_name, "_DE")]][gene_name,])
  
  pos1 = grep("t_adj", colnames(rdata))
  data.frame("model"=sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) x[[1]]),
             "comparison"=sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) {
               xlen = length(x)
               paste(x[c(xlen-1,xlen)], collapse="_")
             }),
             "sex"=sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) {
               xlen = length(x)
               x[[xlen-1]]
             }),
             "group1"=sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) {
               xlen = length(x)
               unlist(strsplit(x[[xlen]], "\\."))[[1]]
             }),
             "group2"=sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) {
               xlen = length(x)
               unlist(strsplit(x[[xlen]], "\\."))[[2]]
             }),
             "domain"= sapply(strsplit(colnames(rdata)[pos1], "_"), function(x) {
               xlen = length(x)
               if(xlen==6) return(x[[4]])
               return(paste0("all",substr(annot_name, 7, 10)))
             }),
             "t_adj.P.Val"=as.numeric(rdata[1,pos1])
  )
}


# format t-test results for ggpubr
annotStandin <- function(logcount_DF, tstat_DF) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  # create standin data with the formatting and attributes expected from ggpubr
  standin = logcount_DF %>% group_by(domain, sex) %>%
    t_test(plot.gene ~ condition)
  
  for(i in 1:nrow(standin)) {
    standin[i,"p.adj"] = filter(tstat_DF, sex==standin$sex[[i]], domain==standin$domain[[i]], 
                                group1==standin$group1[[i]], group2==standin$group2[[i]])$t_adj.P.Val
  }
  standin$p.adj.signif = ifelse(standin$p.adj<.05, "*", "ns")
  
  return(standin)
}

# create empty dataframe if not present in model (needed for whole tissue plots only)
annotEmpty <- function(logcount_DF) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  standin = logcount_DF %>% group_by(domain, sex) %>%
    t_test(plot.gene ~ condition)
  standin$p.adj.signif = "ns"
  
  return(standin)
}

# custom y position adjustments
adjustYposition <- function(logcount_DF, stat_DF) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  #reset the baseline
  adjust.step1 = filter(stat_DF, p.adj.signif!="ns") 
  if(nrow(adjust.step1)==0) return(stat_DF)
  adj1 = group_by(logcount_DF, domain, sex) %>% summarise(max1=max(plot.gene)) %>%
    right_join(stat_DF, by=c("sex","domain")) %>%
    select(domain, sex, max1, y.position)
  stat_DF$y.position = adj1$max1+.5
  
  # now if there are only 2 sig per group, introduce step
  adjust.step2 = filter(stat_DF, p.adj.signif!="ns") %>% 
    group_by(sex, domain) %>% tally() %>%
    filter(n==2)
  if(nrow(adjust.step2)>0) {
    for(i in 1:nrow(adjust.step2)) {
      cond_row = stat_DF$sex==adjust.step2$sex[[i]] & stat_DF$domain==adjust.step2$domain[[i]] &
        stat_DF$p.adj.signif=="*"
      cond_row[cond_row==T] = c(FALSE,TRUE)
      start.position = stat_DF$y.position[cond_row]
      stat_DF[cond_row, "y.position"] = start.position+.8
    }
  }
  
  # now if there are only 3 sig per group, introduce step for both
  ## i actually haven't found an example where this is the case so i haven't been able to trouble shoot this test
  adjust.step3 = filter(stat_DF, p.adj.signif!="ns") %>% 
    group_by(sex, domain) %>% tally() %>%
    filter(n==3)
  if(nrow(adjust.step3)>0) {
    for(i in 1:nrow(adjust.step3)) {
      cond_row1 = stat_DF$sex==adjust.step3$sex[[i]] & stat_DF$domain==adjust.step3$domain[[i]] & 
        stat_DF$group1=="NTC" & stat_DF$group2=="BPD"
      start.position1 = stat_DF$y.position[cond_row1]
      stat_DF[cond_row1, "y.position"] = start.position1+.8
      
      cond_row2 = stat_DF$sex==adjust.step3$sex[[i]] & stat_DF$domain==adjust.step3$domain[[i]] & 
        stat_DF$group1=="MDD" & stat_DF$group2=="BPD"
      start.position2 = stat_DF$y.position[cond_row2]
      stat_DF[cond_row1, "y.position"] = start.position2+1.6
    }
  }
  return(stat_DF)
}


# putting it all together: domain restricted, revised for iSEE
DOMAIN_RESTRICTED <- function(se, gene1, annot_name) {
  # get pseudobulk count info
  log.df = extractLogcounts(gene1, se, annot_name)
  cross.df = crossbarLogcounts(log.df)
  
  if(!inModel(gene1, se, annot_name)) {
    
    p1 <- plotDomainRestricted(log.df, cross.df, annot_name, se, ymax1= ceiling(max(log.df[,5])))
  
  } else {
    #get stats
    t.df = getTstats(gene1, se, annot_name)
    
    stat.df = annotStandin(log.df, t.df) %>%
      add_xy_position(x="domain", dodge=.8, scales="fixed", step.increase=0)
    stat.df = adjustYposition(log.df, stat.df)
    
    
    p1 <- plotDomainRestricted(log.df, cross.df, annot_name, se, ymax1= max(ceiling(c(max(log.df[,5]), max(stat.df$y.position)))))+
      stat_pvalue_manual(stat.df, label="p.adj.signif", hide.ns=T, label.size = 8,
                         color=ifelse(isFsig(gene1, log.df, se, annot_name),"black","grey50"))
  }
  
  return(p1)
}

# putting it all together: whole tissue, revised for iSEE
WHOLE_TISSUE <- function(se, gene1) {
  #function(se, rows, columns) {
  # domain-SP first
  log.df_sp = extractLogcounts(gene1, se, "domain-SP") %>%
    mutate(domain="all-SP")
  cross.df_sp = crossbarLogcounts(log.df_sp)
  # if in model, get stats
  if(!inModel(gene1, se, "domain-SP")) {
    stat.df_sp = annotEmpty(log.df_sp) %>%
      add_xy_position(x="sex", dodge=.8, scales="fixed", step.increase=0)
    stat.df_sp = adjustYposition(log.df_sp, stat.df_sp)
  } else {
    t.df_sp = getTstats(gene1, se, "domain-SP")
    stat.df_sp = annotStandin(log.df_sp, t.df_sp) %>%
      add_xy_position(x="sex", dodge=.8, scales="fixed", step.increase=0)
    stat.df_sp = adjustYposition(log.df_sp, stat.df_sp)
  }
  
  # domain-CT next
  log.df_ct = extractLogcounts(gene1, se, "domain-CT") %>%
    mutate(domain="all-CT")
  cross.df_ct = crossbarLogcounts(log.df_ct)
  # if in model get stats
  if(!inModel(gene1, se, "domain-CT")) {
    stat.df_ct = annotEmpty(log.df_ct) %>%
      add_xy_position(x="sex", dodge=.8, scales="fixed", step.increase=0)
    stat.df_ct = adjustYposition(log.df_ct, stat.df_ct)
  } else {
    t.df_ct = getTstats(gene1, se, "domain-CT")
    stat.df_ct = annotStandin(log.df_ct, t.df_ct) %>%
      add_xy_position(x="sex", dodge=.8, scales="fixed", step.increase=0)
    stat.df_ct = adjustYposition(log.df_ct, stat.df_ct)
  }
  
  log.df = rbind(log.df_sp, log.df_ct) %>% mutate(domain=factor(domain, levels=c("all-SP","all-CT")))
  stat.df = rbind(stat.df_sp, stat.df_ct) %>% mutate(domain=factor(domain, levels=c("all-SP","all-CT")))
  cross.df = rbind(cross.df_sp, cross.df_ct) %>% mutate(domain=factor(domain, levels=c("all-SP","all-CT")))
  
  p1 <- plotWholeTissue(log.df, cross.df, spe_pseudo=se, ymax1= max(ceiling(c(max(log.df[,5]), max(stat.df$y.position), na.rm=T))))
  
  p1+stat_pvalue_manual(filter(stat.df, domain=="all-SP"), label="p.adj.signif", hide.ns=T, label.size = 8,
                        color= ifelse(isFsig(gene1, log.df_sp, se, "domain-SP"),"black","grey50"))+
    stat_pvalue_manual(filter(stat.df, domain=="all-CT"), label="p.adj.signif", hide.ns=T, label.size = 8,
                       color= ifelse(isFsig(gene1, log.df_ct, se, "domain-CT"),"black","grey50"))
  
}

EMPTY <- function(mode) {
  if(mode=="whole-tissue") {
    p1 <- ggplot()+labs(y="logcounts", title="Use drop-down bar to select a gene to plot.")+
      facet_wrap(factor(c("domain-SP","domain-CT"), levels=c("domain-SP","domain-CT")))+
      theme_bw()+theme(text=element_text(size=16))
  }
  if(mode=="domain-restricted") {
    p1 <- ggplot()+labs(y="logcounts", title="Use drop-down bar to select a gene to plot.")+
      facet_wrap(factor(c("Female","Male"), levels=c("Female","Male")))+
      theme_bw()+theme(text=element_text(size=16))
  }
  return(p1)
}
