#  File ergm/R/ergm.geodistn.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
# Code for calculating the geodesic distribution
# 12/03/2004  DH

ergm.geodistn <- function(edgelist, n=max(edgelist), directed=FALSE) {
# edgelist is an mx2 matrix of edges.  n is the number of nodes.
# This function returns a vector of length n, where
#       v[i], i=1, ..., n-1 :  # of pairs of geodesic length i
#       v[n]  : # of pairs of geodesic length infinity.
# Note:  This code does very little error-checking, so don't screw it up
# with illegal vertex numbers (non-positive integers) or an illegal value
# of n.
  
  if(!directed){
   ndyads <- n*(n-1)/2
  }else{
   ndyads <- n*(n-1)
  }
# The C code requires the edgelist to be directed and sorted correctly.
#
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



ergm.geodesicmatrix.edgelist <- function(edgelist, n=max(edgelist), directed=FALSE) {
# edgelist is an mx2 matrix of edges.  n=#nodes.
# This function returns an nxn matrix, where
#      M[i,j] : length of shortest path from vertex i to vertex j
  
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

ergm.nodegeodesics <- function(edgelist, s, n=max(edgelist), directed=FALSE) {
# edgelist is an mx2 matrix of edges.  s=source node. n=#nodes.
# This function returns a vector of length n, where
#      v[i] : length of shortest path from vertex s to vertex i
  
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

ergm.pairgeodesic <- function(edgelist, s, d, n=max(edgelist), directed=FALSE) {
# edgelist is an mx2 matrix of edges.  s=source. d=destination. n=#nodes.
# This function returns the length of the geodesic from s to d.
  
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

