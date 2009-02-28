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

ergm.geodistdist<-function(nw, directed=is.directed(nw)){
 ergm.geodistn(edgelist=as.matrix.network(nw,matrix.type="edgelist"),
               n=nw$gal$n, directed=directed)/(2-is.directed(nw))
}

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



