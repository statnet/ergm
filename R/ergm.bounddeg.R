#  File ergm/R/ergm.bounddeg.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################
# The <ergm.bounddeg> function initializes the list of parameters used
# to bound the degree during the sampling process, and issues warnings
# if the original network doesn't meet the constraints specified by
# 'bounddeg'
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

