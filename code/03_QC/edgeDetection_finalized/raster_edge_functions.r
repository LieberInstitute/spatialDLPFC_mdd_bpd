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

clumpEdges <- function(.xyz) { #.xyz is a dataframe with rownames equal to spot codes with array_row, array_col, then the binary outlier variable
 if(sum(.xyz[,3])==0) return(c())
 t1 = rasterFromXYZ(.xyz)
 rast = as.matrix(t1)
 c1 = clump(t1, directions=8)
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
