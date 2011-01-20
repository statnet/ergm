#===========================================================
# This file contains the 2 following functions for updating
# networks with a proposed set of changes 
#     <network.toggle>
#     <network.accumulate>
#===========================================================



################################################################
# The <network.toggle> and <network.accumulate> functions each
# realize a set of proposed changes on a given network. The
# toggle variant toggles the given edges. The accumulate variant
# adds ?? the given edges
#
# --PARAMETERS--
#   nw      : the network object that will receive the changes
#   nws     : the edges to be toggled or accumulated, as either 
#             a network, an edgelist or a network series; if
#             'timestep' is provided, only a network series is 
#             valid input; if 'timestep' is not provided, only 
#             the edgelist and network are valid input
#   timestep: an indicator of which input types are allowd by
#             'nws'; NULL allows networks and edgelists;
#             non-NULL allows network series; default=NULL
#
# --RETURNED--
#   the original 'nw' with the proposed set of edges toggled or
#   accumulated according to the funtion call
#
#################################################################

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
   nws <- nws$changed[timesteps %in% timestep,2:3,drop=FALSE]
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




#############################################################
# The <network.accumulate> is nearly identical to the
# <network.toggle> function.  please see its function header
#############################################################

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
   nws <- nws$changed[timesteps %in% timestep,2:3,drop=FALSE]
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
