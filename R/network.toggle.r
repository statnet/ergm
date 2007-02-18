network.toggle<-function(nw,diffedges)
{
  if(is.network(diffedges)){
    diffedges <- as.matrix.network(diffedges,"edgelist")
  }else{
    matrix.type <- which.matrix.type(diffedges)
    if(nrow(diffedges)==0){matrix.type <- "edgelist"}
    if(matrix.type!="edgelist"){ 
      stop("network.update requires an edgelist")
    }
  }
  for(i in 1:nrow(diffedges)){  
    nw[diffedges[i,1],diffedges[i,2]] <- 1-nw[diffedges[i,1],diffedges[i,2]] 
  }
  nw
}
