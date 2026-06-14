suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(ggbeeswarm)
	library(rstatix)
	library(ggpubr)
	library(gridExtra)
})

name2id <- function(gene_name, spe_pseudo) {
  gid = rownames(rowData(spe_pseudo))[rowData(spe_pseudo)[,"gene_name"]==gene_name]
  names(gid) = gene_name
  return(gid)
}


# whole tissue dataframe
extractLogcounts <- function(gene_id, spe_pseudo) {
  spe_pseudo[[names(gene_id)]] = logcounts(spe_pseudo)[gene_id,]
  df = as.data.frame(colData(spe_pseudo)[,c("sample_id","condition","sex","domain",names(gene_id))])
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


# create ggplot base
plotFunction <- function(logcount_DF, summary_DF, annot_name, spe_pseudo) {
  save_name = colnames(logcount_DF)[5]
  colnames(logcount_DF)[5] = "plot.gene"
  
  # get limits (just ceiling of max for now, will need to change later with ggpubr)
  ymax1 = ceiling(max(logcount_DF$plot.gene))
  
  # subtitle formating based on sig F test
  subtitle_face = ifelse(isFsig(save_name, logcount_DF, spe_pseudo),"bold","plain")
  # subtitle text to include model name
  if(length(unique(logcount_DF$domain))>1) {
    sub_text="Domain-restricted DE"
  } else {
    sub_text="Whole-tissue DE"
  }

  dx.pal = c("NTC"="#7f7f7f","MDD"="#FB8861","BPD"="#9260b2")

  ggplot(logcount_DF, aes(x=domain))+
    geom_quasirandom(aes(y=plot.gene, color=condition), dodge.width=.8, cex=1)+
    scale_color_manual("DX", values=dx.pal)+
    geom_crossbar(data=summary_DF, aes(y=yavg, ymin=y_min, ymax=y_max, group=condition),
                  color="black", fill="transparent", position=position_dodge(width=.8), width=.6)+
    facet_grid(cols=vars(sex), labeller= as_labeller(c("F"="Female","M"="Male")))+
    ylim(0,ymax1)+
    labs(title=save_name, subtitle=paste(sub_text, getFname(save_name, logcount_DF, spe_pseudo), sep="\n"), 
         y="logcounts", x=annot_name)+
    theme_minimal()+theme(plot.title=element_text(face="italic"), plot.subtitle = element_text(face=subtitle_face))
}


# extract F statistics for labeling
getFname <- function(gene_name, logcount_DF, spe_pseudo) {
  rdata= rowData(spe_pseudo)[name2id(gene_name, spe_pseudo),]
  
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
isFsig <- function(gene_name, logcount_DF, spe_pseudo) {
  rdata= rowData(spe_pseudo)[name2id(gene_name, spe_pseudo),]
  
  if(length(unique(logcount_DF$domain))==1) {
    fpadj = rdata$whole.tissue_F_adj.P.Val
  } else {
    fpadj = rdata$domain.restricted_F_adj.P.Val
  }
  
  return(fpadj<.05)
}

# extract t-test results
getTstats <- function(gene_name, spe_pseudo) {
  rdata= as.data.frame(rowData(spe_pseudo)[name2id(gene_name, spe_pseudo),])
  
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
               return("all")
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
  
  stat_DF = add_xy_position(standin, x="domain", dodge=.8, scales="fixed", step.increase=0)
  stat_DF = adjustYposition(logcount_DF, stat_DF)
  return(stat_DF)
}

# custom y position adjustments
adjustYposition <- function(logcount_DF, stat_DF) {
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


# putting it all together
iSEEplots <- function(spe_pseudo, gene1, ...) {
  # get counts and stats
  log.df = extractLogcounts(name2id(gene1, spe_pseudo), spe_pseudo)
  log.df_all = mutate(log.df, domain="all")
  
  t.df = getTstats(gene1, spe_pseudo)
  
  stat.df = annotStandin(log.df, t.df)
  stat.df_all = annotStandin(log.df_all, t.df)
  
  ymax1 = max(ceiling(c(max(log.df[,5]), max(stat.df$y.position), 
                        max(stat.df_all$y.position))), na.rm=T)
  
  # summarise counts for cross bars 
  cross.df = crossbarLogcounts(log.df)
  cross.df_all = crossbarLogcounts(log.df_all)
  
  # plot for domain-restricted
  p1 <- plotFunction(log.df, cross.df, spe_pseudo, ...)+
    stat_pvalue_manual(stat.df, label="p.adj.signif", hide.ns=T, label.size = 6,
                       color=ifelse(isFsig(gene1, log.df, spe_pseudo),"black","grey50"))+
    ylim(0,ymax1)
  
  # plot for whole-tissue
  p2 <- plotFunction(log.df_all, cross.df_all, spe_pseudo, ...)+
    stat_pvalue_manual(stat.df_all, label="p.adj.signif", hide.ns=T, label.size = 6,
                       color= ifelse(isFsig(gene1, log.df_all, spe_pseudo),"black","grey50"))+
    ylim(0,ymax1)+
    theme(legend.position="none")
  
  
  return(grid.arrange(p2, p1, layout_matrix=matrix(c(1,2,2,2), ncol=4)))
}
