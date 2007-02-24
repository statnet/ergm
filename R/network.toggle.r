network.toggle<-function(nw,nws,time=NULL)
{
  if(is.null(time)){
   if(is.network(nws)){
    nws <- as.matrix.network(nws,"edgelist")
   }else{
    matrix.type <- which.matrix.type(nws)
    if(nrow(nws)==0){matrix.type <- "edgelist"}
    if(matrix.type!="edgelist"){ 
      stop("network.toggle requires an edgelist or a network series")
    }
   }
  }else{
   if(class(nws)!="network.series"){ 
      stop("network.toggle requires an edgelist or a network series")
   }
   times <- nws$changed[,1]
   nws <- nws$changed[times==time,2:3]
  }
  for(i in 1:nrow(nws)){  
    nw[nws[i,1],nws[i,2]] <- 1-nw[nws[i,1],nws[i,2]] 
  }
  nw
}
