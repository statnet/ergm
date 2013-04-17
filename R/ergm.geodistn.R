#  File R/ergm.geodistn.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#======================================================================================
# This file contains the following 6 functions for computing various geodesic measures
#        <ergm.geodistdist>         <ergm.geodesicmatrix.edgelist>
#        <ergm.geodistn>            <ergm.nodegeodesics>
#        <ergm.geodesicmatrix>      <ergm.pairgeodesic>
#======================================================================================



ergm.geodistdist<-function(nw, directed=is.directed(nw)){
 ergm.geodistn(edgelist=as.edgelist(nw),
               n=network.size(nw), directed=directed)/(2-is.directed(nw))
}


###############################################################################
# The <ergm.geodistn> function calculates and returns the geodesic  distance
# distribution for a given network via <full_geodesic_distribution.C>
# Note:  This code does very little error-checking, so don't screw it up
# with illegal vertex numbers (non-positive integers) or an illegal value
# of n.
#
# --PARAMETERS--
#   edgelist:  the edgelist an mx2 matrix
#   n       :  the number of nodes in the network; default=max(edgelist)
#   directed:  whether the edgelist represents a directed network (T or F);
#              default=FALSE
#
#
# --RETURNED--
#   ans: an n-length vector where
#      ans[i], i=1, ..., n-1 is the number of pairs of geodesic length i
#      ans[n] is the number of pairs of geodesic length infinity.
#
################################################################################

ergm.geodistn <- function(edgelist, n=max(edgelist), directed=FALSE) {
  if(!directed){
   ndyads <- n*(n-1)/2
  }else{
   ndyads <- n*(n-1)
  }
# The C code requires the edgelist to be directed and sorted correctly.
  if(!is.matrix(edgelist) || nrow(edgelist)==0){
   return(rep(c(0,ndyads),c(n-1,1)))
  }
  if(nrow(edgelist)>1){
   edgelist<-edgelist[edgelist[,1]!=edgelist[,2],] # get rid of self-edges
   if (!directed) 
    edgelist<-rbind(edgelist,edgelist[,2:1])
   edgelist<-unique(edgelist)
   edgelist<-edgelist[order(edgelist[,1],edgelist[,2]),]
  }else{
   if(edgelist[1]==edgelist[2]){return(rep(c(0,ndyads),c(n-1,1)))}
   return(rep(c(1,0,ndyads-1),c(1,n-2,1)))
  }

# Next, we need to set up the nodelist vector:  Because of C's numbering
# convention, we want nodelist[1]=0 and in general, nodelist[i]=2*r(i)-2,
# where r(i) is the first row in edgelist containing from node i.  (If
# there are no edges from node i, just set nodelist[i]=0.)
  nodelist<-match(1:n,edgelist[,1],nomatch=1)-1
  
# Now everything is ready.  Call the C code.
  ans<-.C("full_geodesic_distribution", as.integer(t(edgelist)),
    as.integer(n), as.integer(nodelist), as.integer(dim(edgelist)[1]),
    colors=integer(n), distances=integer(n), queue=integer(n),
    distribution=integer(n), PACKAGE='ergm') $ distribution
  names(ans)<-c(1:(n-1),"Inf") # length n really means no path exists
  ans
}




ergm.geodesicmatrix <- function(nw, directed=is.directed(nw)){
 ergm.geodesicmatrix.edgelist(edgelist=as.edgelist(nw),
               n=network.size(nw), directed=directed)
}



##################################################################################
# The <ergm.geodesicmatrix.edgelist> function calculates and returns a matrix of
# the shorted path lengths between all node for a given network, nw, via
# <geodesic_matrix.C>
#
# --PARAMETERS--
#   edgelist:  the edgelist of the network  
#   n       :  the number of nodes in the network; default=max(edgelist) 
#   directed:  whether the edgelist represents a directed network (T or F);
#              default=FALSE
#
# --RETURNED--
#   an n x n matrix, whose i,j entry is the shortest path length between nodes
#   i and j
#
###################################################################################

ergm.geodesicmatrix.edgelist <- function(edgelist, n=max(edgelist), directed=FALSE) {
# This function starts off just like ergm.geodistn:
  edgelist<-edgelist[edgelist[,1]!=edgelist[,2],] # get rid of self-edges
  if (!directed) 
    edgelist<-rbind(edgelist,edgelist[,2:1])
  edgelist<-unique(edgelist)
  edgelist<-edgelist[order(edgelist[,1],edgelist[,2]),]
  nodelist<-match(1:n,edgelist[,1],nomatch=1)-1
  
# Now everything is ready.  Call the C code.
  ans<-.C("geodesic_matrix", as.integer(t(edgelist)), as.integer(n),
    as.integer(nodelist), as.integer(dim(edgelist)[1]), colors=integer(n),
    gmat=integer(n*n), queue=integer(n), PACKAGE='ergm') $ gmat
  ans[ans==n]<-Inf # length n really means no path exists
  ans=matrix(ans,n,n,byrow=TRUE) # byrow=TRUE is only important when directed==TRUE
  ans
}






##################################################################################
# The <ergm.nodegeodesics> function calculates and returns a vector of the
# shortest path lengths between a given node s and all others of a network via
# <node_geodesics.C>
#
# --PARAMETERS--
#   edgelist:  the edgelist of the network
#   s       :  the node from which to calculate shortest paths; 1<s<n
#   n       :  the number of nodes in the network; default=max(edgelist)
#   directed:  whether the edgelist represents a directed network (T or F);
#              default=FALSE
#
# --RETURNED--
#   ans: a vector of length n whose ith entry is the length of shortest path from
#        vertex s to vertex i
#
###################################################################################

ergm.nodegeodesics <- function(edgelist, s, n=max(edgelist), directed=FALSE) {
  
# This function starts off just like ergm.geodistn:
  edgelist<-edgelist[edgelist[,1]!=edgelist[,2],] # get rid of self-edges
  if (!directed) 
    edgelist<-rbind(edgelist,edgelist[,2:1])
  edgelist<-unique(edgelist)
  edgelist<-edgelist[order(edgelist[,1],edgelist[,2]),]
  nodelist<-match(1:n,edgelist[,1],nomatch=1)-1
  
# Now everything is ready.  Call the C code.
  ans<-.C("node_geodesics", as.integer(t(edgelist)), as.integer(n),
    as.integer(nodelist), as.integer(dim(edgelist)[1]), colors=integer(n),
    distances=integer(n), queue=integer(n), as.integer(s), PACKAGE='ergm') $ distances
  ans[ans==n]<-Inf # length n really means no path exists
  ans
}





##################################################################################
# The <ergm.pairgeodesic> function calculates and returns the geodesic distance
# between a pair of nodes, s and d, in a network, via <pair_geodesic.C>
#
# --PARAMETERS--
#   edgelist:  the edgelist for the network
#   s       :  the source node 
#   d       :  the destination node
#   n       :  the number of nodes in the network; default=max(edgelist)
#   directed:  whether the edgelist represents a directed network (T or F);
#              default=FALSE
#
# --RETURNED--
#   ans: a vector of length n whose ith entry is the length of shortest path from
#        vertex s to vertex i
#
###################################################################################

ergm.pairgeodesic <- function(edgelist, s, d, n=max(edgelist), directed=FALSE) {
# This function starts off just like ergm.geodistn:
  edgelist<-edgelist[edgelist[,1]!=edgelist[,2],] # get rid of self-edges
  if (!directed) 
    edgelist<-rbind(edgelist,edgelist[,2:1])
  edgelist<-unique(edgelist)
  edgelist<-edgelist[order(edgelist[,1],edgelist[,2]),]
  nodelist<-match(1:n,edgelist[,1],nomatch=1)-1
  
# Now everything is ready.  Call the C code.
  ans<-.C("pair_geodesic", as.integer(t(edgelist)), as.integer(n),
    as.integer(nodelist), as.integer(dim(edgelist)[1]), colors=integer(n),
    distances=integer(n), queue=integer(n), as.integer(s),
    as.integer(d), PACKAGE='ergm') $ distances[d]
  if (ans==n) ans<-Inf # length n really means no path exists
  ans
}

