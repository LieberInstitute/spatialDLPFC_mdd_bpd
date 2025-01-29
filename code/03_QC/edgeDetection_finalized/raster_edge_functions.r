my_fill <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  x[is.na(x)] = 1
  if(sum(x[-5])==8) return(1)
  else return(0)
}
my_fill_star <- function(x) {
  if(is.na(x[5])) return(NA_integer_)
  if(x[5]==1) return(1)
  if(sum(x[-5], na.rm=T)==4) return(1)
  else return(0)
}

my_outline <- function(x) {
  if(is.na(x[13])) return(NA_integer_)
  if(x[13]==1) return(1)
  x[is.na(x)] = 1
  sum.edges = sum(x[1:5])+sum(x[21:25])+sum(x[c(6,11,16)])+sum(x[c(10,15,20)])
  if(sum.edges==16) return(1)
  else(return(0))
}

lookupKey <- function(.xy, results_df) { #.xy is a dataframe with rownames equal to spot codes with array_row, array_col
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  t1 = rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = sapply(1:nrow(results_df), function(x) {
    v1 = as.numeric(results_df[x,])
    rast[v1[1],v1[2]]
  })
  edge_spots = rownames(key1)[key1$numeric_key %in% num_key] 
  return(edge_spots)
}

clumpEdges <- function(.xyz, shifted=FALSE) { #.xyz is a dataframe with rownames equal to spot codes with array_row, array_col, then the binary outlier variable
 if(sum(.xyz[,3])==0) return(c())
 if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
 }
 t1 = rasterFromXYZ(.xyz)
 #fill any lone empty spots in a sea out outliers
 t2 <- focal(t1, w=matrix(1,3,3), fun=my_fill, pad=T, padValue=1)
 rast = as.matrix(t2)
 c1 = clump(t2, directions=8)
 clumps = as.matrix(c1)
 
 edgeClumps = c()
 north = sum(!is.na(clumps[1,]))/sum(!is.na(rast[1,]))
 if(north>=.75) edgeClumps = c(edgeClumps, unique(clumps[1,]))
 east = sum(!is.na(clumps[,ncol(clumps)]))/sum(!is.na(rast[,ncol(rast)]))
 if(east>=.75) edgeClumps = c(edgeClumps, unique(clumps[,ncol(clumps)]))
 south = sum(!is.na(clumps[nrow(clumps),]))/sum(!is.na(rast[nrow(rast),]))
 if(south>=.75) edgeClumps = c(edgeClumps, unique(clumps[nrow(clumps),]))
 west = sum(!is.na(clumps[,1]))/sum(!is.na(rast[,1]))
 if(west>=.75) edgeClumps = c(edgeClumps, unique(clumps[,1]))

 edgeClumps = edgeClumps[!is.na(edgeClumps)]
 if(length(edgeClumps)==0) return(c())
 
 res <- vector("list", length(edgeClumps))
 names(res) <- as.character(edgeClumps)
 for (i in edgeClumps){
   res[[as.character(i)]] <- as.data.frame(which(clumps == i, arr.ind = TRUE))
 }
 res.df = do.call(rbind, res) 
 
 return(lookupKey(.xyz[,1:2], res.df))
}

lookupKeySize <- function(.xy, results_df) { #.xyz is a dataframe with rownames equal to spot codes with array_row, array_col and results_df is a sample specific df with 3 columns for row, col, and size of clump
  key1 = cbind(.xy, "numeric_key"=seq(nrow(.xy)))
  #key1 = cbind(test[,1:2], "numeric_key"=seq(nrow(test)))
  t1 = rasterFromXYZ(key1)
  rast = as.matrix(t1)
  
  num_key = do.call(rbind, lapply(1:nrow(results_df), function(x) {
    v1 = as.numeric(results_df[x,])
    cbind.data.frame("numeric_key"=rast[v1[1],v1[2]], "size"=v1[3])
  }))
  matching = match(num_key$numeric_key, key1$numeric_key)
  return(cbind.data.frame("spotcode"=rownames(key1)[matching], "size"=num_key$size))
}

problemAreas <- function(.xyz, shifted=FALSE) { #.xyz is a dataframe with rownames equal to spot codes with array_row, array_col, then the binary outlier variable
  if(sum(.xyz[,3])==0) return(c())
  if(shifted==TRUE) {
    odds = seq(1,max(.xyz[,"array_col"]), by=2)
    .xyz[.xyz[,"array_col"] %in% odds, "array_col"] = .xyz[.xyz[,"array_col"] %in% odds, "array_col"]-1
  }
  t1 = rasterFromXYZ(.xyz)
  t2 <- focal(t1, w=matrix(1,3,3), fun=my_fill, pad=T, padValue=1)
  t3 <- focal(t2, w=matrix(1,5,5), fun=my_outline, pad=T, padValue=1)
  star.m = matrix(rep(c(0,1), length.out=9), 3, 3)
  star.m[5] = 1
  t3_s <- focal(t3, star.m, fun=my_fill_star, pad=T, padValue=1)
  rev_t3 = t3_s
  rev_t3[t3_s==0] = 1
  rev_t3[t3_s==1] = 0
  rev_c3 = clump(rev_t3)
  tbl = table(as.matrix(rev_c3))
  not_clump = as.numeric(names(tbl)[tbl==max(tbl)])
  t4 = t3_s
  t4[rev_c3!=not_clump] = 1
  #t1 = rasterFromXYZ(m1)
  rast = as.matrix(t4)
  c1 = clump(t4, directions=8) #plot c1 without shifting array col doesn't produce continuous tissue thread
  clumps = as.matrix(c1)
  
  tot <- max(clumps, na.rm=TRUE)
  res <- vector("list",tot)
  for (i in 1:tot){
    res[i] <- list(which(clumps == i, arr.ind = TRUE))
    res[[i]] <- cbind(res[[i]],"size"=nrow(res[[i]]))
  }

  res.df = do.call(rbind, res) 
  
  return(lookupKeySize(.xyz[,1:2], res.df))
}
