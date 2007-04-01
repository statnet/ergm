network.toggle<-function(nw,nws,timestep=NULL)
{
  if(is.null(timestep)){
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
   timesteps <- nws$changed[,1,drop=FALSE]
   nws <- nws$changed[timestep %in% timesteps,2:3,drop=FALSE]
  }
  if(nrow(nws) > 0){
#  for(i in 1:nrow(nws)){  
#   nw[nws[i,1],nws[i,2]] <- 1-nw[nws[i,1],nws[i,2]] 
#  }
   toggle.dyads(nw,tail=nws[,1],head=nws[,2])
  }else{
   nw
  }
}
network.accumulate<-function(nw,nws,timestep=NULL)
{
  if(is.null(timestep)){
   if(is.network(nws)){
    nws <- as.matrix.network(nws,"edgelist")
   }else{
    matrix.type <- which.matrix.type(nws)
    if(nrow(nws)==0){matrix.type <- "edgelist"}
    if(matrix.type!="edgelist"){ 
      stop("network.accumulate requires an edgelist or a network series")
    }
   }
  }else{
   if(class(nws)!="network.series"){ 
      stop("network.accumulate requires an edgelist or a network series")
   }
   timesteps <- nws$changed[,1,drop=FALSE]
   nws <- nws$changed[timesteps==timestep,2:3,drop=FALSE]
  }
  if(nrow(nws) > 0){
#  for(i in 1:nrow(nws)){  
#   nw[nws[i,1],nws[i,2]] <- 1-nw[nws[i,1],nws[i,2]] 
#  }
   accumulate.edges(nw,tail=nws[,1],head=nws[,2])
  }else{
   nw
  }
}
