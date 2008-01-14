ergm.boundDeg <- function(boundDeg,nw){    
  #  Resolve conditioning in ERGM call, as expressed in the
  #  argument boundDeg (a list, with item names as seen below)
  nnodes=network.size(nw)
  if(is.null(boundDeg) ||
     all(sapply(boundDeg,function(x){length(x)=1 && x==0}))
  ) {
    attribs <- 0
    maxout <- 0
    maxin <- 0
    minout <- 0
    minin <- 0
  } else {
    attribs <- boundDeg$attribs
    maxout <- boundDeg$maxout
    maxin <- boundDeg$maxin
    minout <- boundDeg$minout
    minin <- boundDeg$minin
    if (is.null(attribs) || attribs==0){ 
      if(any(!is.null(c(minin,minout,maxout,maxin)))){ 
        attribs <- matrix(1,ncol=1,nrow=nnodes)
        # Get degree for each node
        el <- as.matrix.network.edgelist(nw)
        if (!is.directed(nw)) {
          outdeg <- tabulate(as.vector(el), nbins=nnodes)
        } else {
          outdeg <- tabulate(el[,1], nbins=nnodes)
          indeg <- tabulate(el[,2], nbins=nnodes)
        }
      }else{
        attribs <- 0
      }
    }
    if(is.null(minin )) minin <- matrix(0,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(minout)) minout <- matrix(0,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(maxin )) maxin <- matrix(nnodes-1,ncol=ncol(attribs),nrow=nnodes)
    if(is.null(maxout)) maxout <- matrix(nnodes-1,ncol=ncol(attribs),nrow=nnodes)
    if(length(minin )==1) minin  <- matrix(minin ,ncol=ncol(attribs),nrow=nnodes)
    if(length(minout)==1) minout <- matrix(minout,ncol=ncol(attribs),nrow=nnodes)
    if(length(maxin )==1) maxin  <- matrix(maxin ,ncol=ncol(attribs),nrow=nnodes)
    if(length(maxout)==1) maxout <- matrix(maxout,ncol=ncol(attribs),nrow=nnodes)
    minin[is.na( minin)| minin<0] <- 0
    minout[is.na(minout)|minout<0] <- 0
    maxin[is.na( maxin)] <- nnodes-1
    maxout[is.na(maxout)] <- nnodes-1
    if (any(outdeg>maxout | outdeg<minout) || (is.directed(nw) && any(indeg>maxin | indeg<minin))) {
      cat("Warning!  Initial network does not satisfy degree constraints.\n",
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

