#  File R/ergm.geodistn.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
#======================================================================================
# This file contains the following 6 functions for computing various geodesic measures
#        <ergm.geodistdist>         <ergm.geodesicmatrix.edgelist>
#        <ergm.geodistn>            <ergm.nodegeodesics>
#        <ergm.geodesicmatrix>      <ergm.pairgeodesic>
#======================================================================================





#' Calculate geodesic distance distribution for a network or edgelist
#' 
#' \code{ergm.geodistdist} calculates geodesic distance distribution for a
#' given \code{\link{network}} and returns it as a vector.
#' 
#' \code{ergm.geodistdist} is a network wrapper for \code{ergm.geodistn}, which
#' calculates and returns the geodesic distance distribution for a given
#' network via full_geodesic_distribution.C
#' 
#' @param nw \code{\link{network}} object over which distances should be
#' calculated
#' @param directed logical, should the network be treated as directed
#' @param edgelist an edgelist representation of a network as an mx2 matrix
#' @param n integer, size of the network
#' @return a vector \code{ans} with length equal to the size of the network
#' where \itemize{
#' \item `ans[i], i=1, ..., n-1` is the number of pairs of
#' geodesic length `i`
#' \item `ans[n]` is the number of pairs of geodesic length
#' infinity.  }
#' @seealso See also the sna package \code{\link[sna]{geodist}} function
#' @examples
#' 
#' data(faux.mesa.high)
#' ergm.geodistdist(faux.mesa.high)
#'
#' @keywords internal
#' @export ergm.geodistdist
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

#' @rdname ergm.geodistdist
#' @description \code{ergm.geodistn} calculates geodesic deistance
#'   distribution based on an input edgelist, and has very little
#'   error checking so should not normally be called by users. The C
#'   code requires the edgelist to be directed and sorted correctly.
#' 
#' @export
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


