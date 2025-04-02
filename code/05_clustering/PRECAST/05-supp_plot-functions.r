reducedDimPlots <- function(object, .dimred) {
  if(.dimred=="PRECAST_1663") {
    object$precast_cluster_f = object$precast_k9_1663_f
    fill_color = precast.colorList[["n1663_k9"]][["colors"]]
    names(fill_color) = precast.colorList[["n1663_k9"]][["annotation"]]
    ptitle= "n=1663"
  }
  if(.dimred=="PRECAST_1079") {
    object$precast_cluster_f = object$precast_k9_1079_f
    fill_color = precast.colorList[["n1079_k9"]][["colors"]]
    names(fill_color) = precast.colorList[["n1079_k9"]][["annotation"]] 
    ptitle= "n=1079"
  }
  
  test_variables = c("sample_id","slide","sum_umi","precast_cluster_f")#, 
                     #"sex","condition","PMI","RIN")
  percVar.df = tibble::rownames_to_column(as.data.frame(getExplanatoryPCs(object, dimred=.dimred, n_dimred = 15, 
                                                                          variables=colData(object)[,test_variables])), 
                                          var="component") %>%
    tidyr::pivot_longer(all_of(test_variables), names_to="variable", values_to="percVariance")
  ## which experimental or technical variables are in the top 3 of variance for any component
  ## group_by(percVar.df, component) %>% slice_max(n=3, percVariance) %>% group_by(variable) %>% tally()
  ##  variable            n
  ##1 brnum              15
  ##2 precast_k9_1663    15
  ##3 slide               8
  ##4 sum_umi             7
  percVar.df$component = as.numeric(sapply(strsplit(percVar.df$component, "_"), function(x) x[[2]]))
  
  #spaghetti plot of % var per component
  p1 <- ggplot(mutate(percVar.df, variable=factor(variable, levels=c("precast_cluster_f","sample_id","slide","sum_umi"),
                                        labels=c("cluster","sample","slide","sum UMI"))), 
               aes(x=component, y=percVariance, color=variable))+
    geom_point()+geom_line(aes(group=variable))+
    ylim(0,100)+scale_x_continuous(breaks=c(1,5,10,15))+
    scale_color_manual(values=c("black","#66c2a5","#fc8d62","#8da0cb"))+
    labs(title=paste0("PRECAST embeddings (", ptitle ,"): top variables explaining variance"), 
         x="component", y="% of var explained", color="")+
    theme_minimal()+theme(text=element_text(size=10))
  
  #silhouette plot
  set.seed(123)
  sil.results <- as.data.frame(approxSilhouette(reducedDim(object, .dimred), clusters=object$precast_cluster_f))
  sil.results$closest <- factor(ifelse(sil.results$width > 0, as.character(sil.results$cluster), as.character(sil.results$other)))
  sil.results$closest <- factor(sil.results$closest, levels=levels(object$precast_cluster_f))
  
  p2 <- ggplot(sil.results, aes(x=cluster, y=width, colour=closest))+
    ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=fill_color)+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    labs(x="assigned PRECAST cluster", title=paste0("PRECAST embeddings (", ptitle, "): silhouette"))+
    theme_minimal()+theme(text=element_text(size=10))
  
  #sample embeddings per cluster
  med.df = as.data.frame(cbind.data.frame(colData(object)[,c("sample_id","precast_cluster_f")], reducedDim(object, .dimred))) %>%
    group_by(sample_id, precast_cluster_f) %>% summarise_all(median) %>% 
    tidyr::pivot_longer(all_of(paste("PRECAST",1:15, sep="_")), names_to="component", values_to="median_loading") %>%
    mutate(component=factor(component, levels=paste("PRECAST",1:15, sep="_")))
  
  p3 <- ggplot(med.df, aes(x=precast_cluster_f, y=median_loading, color=precast_cluster_f))+
    ggbeeswarm::geom_quasirandom(size=.1)+scale_color_manual(values=fill_color)+
    facet_wrap(vars(component), ncol=4)+labs(y="per-sample median spot loading")+
    #scale_x_discrete(labels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM","low\nUMI"))+
    theme_minimal()+theme(text=element_text(size=8), panel.border=element_rect(fill=NA, color="black", linewidth=.3),
                          axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1),
                          legend.position="none", plot.margin=unit(c(4,8,4,4),"pt"))
  
  #return
  return(list(p1, p2, p3))
}

locationSpotPlots <- function(object_small, .dimred) {
  #object_small = subset spe
  if(.dimred=="PRECAST_1663") {
    object_small$precast_cluster_f = object_small$precast_k9_1663_f
    fill_color = precast.colorList[["n1663_k9"]][["colors"]]
    names(fill_color) = precast.colorList[["n1663_k9"]][["annotation"]]
    ptitle= "n=1663"
  }
  if(.dimred=="PRECAST_1079") {
    object_small$precast_cluster_f = object_small$precast_k9_1079_f
    fill_color = precast.colorList[["n1079_k9"]][["colors"]]
    names(fill_color) = precast.colorList[["n1079_k9"]][["annotation"]]
    ptitle= "n=1079"
  }
  
  sampleList= unique(object_small$sample_id)
  names(sampleList) = sampleList
  sampleList <- lapply(sampleList, function(x) object_small[,object_small$sample_id==x])

  numpal = rep(NA,9)
  names(numpal) = levels(object_small$precast_cluster_f)
  numpal[["GABA"]] = "black"

  plist <- lapply(sampleList, function(y) {
    p <- plotVisium(y, spots=F, image=T) |> add_ground("precast_cluster_f", stroke=.5)
    p+scale_color_manual(values=numpal, na.value = "transparent")+
      theme(legend.position="none", plot.title=element_text(size=8))
  })
  
  return(plist)
}
