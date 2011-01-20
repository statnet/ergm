####################################################################
# The <network.dyadcount> function counts the number of dyads in a
# given network; the count may include or exclude missing dyads
#
# --PARAMETERS--
#   x      :  a network
#   na.omit:  whether to exclude missing dyads from the returned 
#             count; default=TRUE
#
# --RETURNED--
#   dyads:  the count of non-missing dyads in 'x' if 'na.omit'=TRUE
#           the count of all dyads in 'x' if 'na.omit'=FALSE
#
####################################################################

network.dyadcount<-function(x,na.omit=TRUE){
  if(!is.network(x))
    stop("network.dyadcount requires an argument of class network.")

  nodes <- network.size(x)
  if(is.directed(x)){
     dyads <- nodes * (nodes-1)
  }else{
   if(is.bipartite(x)){
    nactor <- get.network.attribute(x,"bipartite")
    nevent <- nodes - nactor
    dyads <- nactor * nevent
   }else{
    dyads <- nodes * (nodes-1)/2
   }
  }
  if(na.omit){
#
#  Adjust for missing
#
   design <- get.network.attribute(x,"design")
   if(!is.null(design)){
    dyads <- dyads - network.edgecount(design)
   }else{
    design <- get.network.attribute(x,"mClist.design")
    if(!is.null(design)){
     dyads <- dyads - design$nedges
    }else{
     dyads <- dyads - network.naedgecount(x)
    }
   }
  }
  dyads
}
