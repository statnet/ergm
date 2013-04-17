#  File R/ergm.bounddeg.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
########################################################################
# The <ergm.bounddeg> function initializes the list of parameters used
# to bound the degree during the sampling process, and issues warnings
# if the original network doesn't meet the constraints specified by
# 'bounddeg'
# 
#
# --PARAMETERS--
#   bounddeg: a list of parameters which may contain the following for
#             a network of size n nodes:
#      attribs: an nxp matrix, where entry ij is TRUE if node i has
#               attribute j, and FALSE otherwise; default=an nx1 matrix
#               of 1's
#      maxout : an nxp matrix, where entry ij is the maximum number of
#               out degrees for node i to nodes with attribute j;
#               default=an nxp matrix of the value (n-1)
#      maxin  : defined similarly to maxout, but ignored for undirected
#               networks; default=an nxp matrix of the value (n-1)
#      minout : defined similarly to maxout; default=an nxp matrix of 0's
#      minin  : defined similarly to maxout, but ignored for undirected
#               networks; default=an nxp matrix of 0's
#   nw: the orginal network specified to <ergm> in 'formula'
#
# --RETURNED--
#   a list of parameters used to bound degree during sampling
#      condAllDegExact: always FALSE
#      attribs        : as defined above
#      maxout         : as defined above
#      maxin          : as defined above
#      minout         : as defined above
#      minin          : as defined above
#   
########################################################################

ergm.bounddeg <- function(bounddeg,nw){    
  nnodes=network.size(nw)
  if(is.null(bounddeg) ||
     all(sapply(bounddeg,function(x){length(x)==1 && x==0}))) {
    attribs <- NULL
    maxout <- NULL
    maxin <- NULL
    minout <- NULL
    minin <- NULL
  } else {
    attribs <- bounddeg$attribs
    maxout <- bounddeg$maxout
    maxin <- bounddeg$maxin
    minout <- bounddeg$minout
    minin <- bounddeg$minin
    if (is.null(attribs) || all(attribs==0)){ 
      if(any(!is.null(c(minin,minout,maxout,maxin)))){ 
        attribs <- matrix(1,ncol=1,nrow=nnodes)
        # Get degree for each node
        el <- as.edgelist(nw)
        if (!is.directed(nw)) {
          outdeg <- tabulate(as.vector(el), nbins=nnodes)
        } else {
          outdeg <- tabulate(el[,1], nbins=nnodes)
          indeg <- tabulate(el[,2], nbins=nnodes)
        }
      }else{
        attribs <- 0
      }
    }else{   #Tabulate degrees by attributes
      el <- as.edgelist(nw)
      if(!is.directed(nw)){
        el<-rbind(el[,1:2],el[,2:1])    #Need edges going both directions
        outdeg<-apply(attribs,2,function(z){tabulate(el[z[el[,2]],1], nbins=nnodes)})
      }else{
        outdeg<-apply(attribs,2,function(z){tabulate(el[z[el[,2]],1], nbins=nnodes)})
        indeg<-apply(attribs,2,function(z){tabulate(el[z[el[,1]],2], nbins=nnodes)})
      }
    }
    if(is.null(minin )) minin <- matrix(0,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(minout)) minout <- matrix(0,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(maxin ) || maxin==0) maxin <- matrix(nnodes-1,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(maxout) || maxout==0) maxout <- matrix(nnodes-1,ncol=ncol(attribs),nrow=nnodes)
    if(length(minin )==1) minin  <- matrix(minin ,ncol=ncol(attribs),nrow=nnodes)
    if(length(minout)==1) minout <- matrix(minout,ncol=ncol(attribs),nrow=nnodes)
    if(length(maxin )==1) maxin  <- matrix(maxin ,ncol=ncol(attribs),nrow=nnodes)
    if(length(maxout)==1) maxout <- matrix(maxout,ncol=ncol(attribs),nrow=nnodes)
    minin[is.na( minin)| minin<0] <- 0
    minout[is.na(minout)|minout<0] <- 0
    maxin[is.na( maxin)] <- nnodes-1
    maxout[is.na(maxout)] <- nnodes-1
    if (any(outdeg>maxout | outdeg<minout) || (is.directed(nw) && any(indeg>maxin | indeg<minin))) {
      cat("Warning:  Initial network does not satisfy degree constraints.\n",
          "Proceeding anyway, but final network may not satisfy constraints.\n")
    }
    attribs[is.na(attribs)] <- 0
  }
  list(condAllDegExact=FALSE,
       attribs=attribs,
       maxout=maxout,
       maxin=maxin,
       minout=minout,
       minin=minin)
}

