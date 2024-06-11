findEdges <- function(coldata, coord.dim=c("array_col","array_row")) {
  if(coord.dim=="array_col") {
    d1 = filter(as.data.frame(coldata), keep_spots=="image perimeter") %>% group_by(array_col) %>% tally() %>% filter(n>10)
    d1 = filter(as.data.frame(coldata), keep_spots=="image perimeter", array_col %in% d1$array_col) %>% group_by(array_col) %>% summarise(n.umi.outlier=sum(umi_3MAD.outlier_slide), n.genes.outlier=sum(genes_3MAD.outlier_slide), n.total=n(), perc.umi=n.umi.outlier/n.total, perc.genes=n.genes.outlier/n.total) 
  }
  if(coord.dim=="array_row") {
    d1 = filter(as.data.frame(coldata), keep_spots=="image perimeter") %>% group_by(array_row) %>% tally() %>% filter(n>10)
    d1 = filter(as.data.frame(coldata), keep_spots=="image perimeter", array_row %in% d1$array_row) %>% group_by(array_row) %>% summarise(n.umi.outlier=sum(umi_3MAD.outlier_slide), n.genes.outlier=sum(genes_3MAD.outlier_slide), n.total=n(), perc.umi=n.umi.outlier/n.total, perc.genes=n.genes.outlier/n.total)
  }
  if(max(d1$perc.umi) >.8 | max(d1$perc.genes) >.8) {
  return(as.numeric(filter(d1, perc.umi>.8 | perc.genes>.8)[,1]))
  }
  else{return(NA)}
}

idEdge <- function(coldata, array.edge) {
  #subset to only outliers
  c1 = filter(as.data.frame(coldata),umi_3MAD.outlier_slide==TRUE | genes_3MAD.outlier_slide==TRUE)
  #pull coord matrix
  m1 = as.matrix(c1[,c("array_row","array_col")])
  rownames(m1) = c1$key
  #find neighbors within 1 that are also outliers (b/c m1 is only outliers)
  test2 = dnearneigh(m1, d1=1, d2=2.1, bounds=c("GE","LE"))#, row.names=c1$key)
  names(test2) = rownames(m1)
  #make an adj matrix out of these neighbors
  a.mtx = nb2mat(test2, zero.policy=TRUE)
  colnames(a.mtx) = names(test2)
  rownames(a.mtx) = names(test2)
  #set up for iterative approach
  new.seed = rownames(m1)[m1[,names(array.edge)]%in%array.edge] #starting row/col
  exclude.spots = new.seed #names in starting row/col
  #iterate
  while(length(new.seed)>1) {
    r1 = colSums(a.mtx[new.seed,])>0
    new.seed = setdiff(colnames(a.mtx)[r1], exclude.spots)
    exclude.spots = union(new.seed, exclude.spots)
  }
  exclude.spots
}
