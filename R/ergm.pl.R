ergm.pl<-function(Clist, Clist.miss=NULL, m, theta.offset=NULL,
                    maxMPLEsamplesize=1e+6,
                    maxNumDyadTypes=1e+6,
                    verbose=FALSE,
                    compressflag=TRUE) {
  offset <- rep(0,Clist$ndyads)
  bip <- Clist$bipartite
  n <- Clist$n
  maxNumDyadTypes <- min(maxNumDyadTypes,
                         ifelse(bip>0, bip*(n-bip), 
                                ifelse(Clist$dir, n*(n-1), n*(n-1)/2)))
  # May have to think harder about what maxNumDyadTypes should be if we 
  # implement a hash-table approach to compression.  
  z <- .C("MPLE_wrapper",
          as.integer(Clist$heads),    as.integer(Clist$tails),
          as.integer(Clist$nedges), as.integer(Clist$maxpossibleedges),
          as.integer(n), 
          as.integer(Clist$dir),     as.integer(bip),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring), as.character(Clist$snamestring),
          as.double(Clist$inputs),
          y = integer(maxNumDyadTypes),
          x = double(maxNumDyadTypes*Clist$nparam),
          weightsvector = integer(maxNumDyadTypes),
          as.double(offset), compressedOffset=double(maxNumDyadTypes),
          as.integer(maxNumDyadTypes),
          as.integer(maxMPLEsamplesize),
          as.integer(compressflag),
          PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  if (verbose) {
    if (compressflag) cat("Compressed ")
    cat(paste("MPLE covariate matrix has", sum(uvals), "rows.\n"))
  }
  zy <- z$y[uvals]
  wend <- z$weightsvector[uvals]
  xmat <- matrix(z$x, ncol=Clist$nparam, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  rm(z,uvals)
  #
  # Adjust for the offset
  #
  if(any(m$etamap$offsettheta)){
    if(is.null(theta.offset)){
      theta.offset <- rep(0, length=Clist$nparam)
      names(theta.offset) <- m$coef.names
      theta.offset[m$etamap$offsettheta] <- -Inf
    }
    foffset <- xmat[,!m$etamap$offsettheta,drop=FALSE]%*%theta.offset[!m$etamap$offsettheta]
    shouldoffset <- apply(abs(xmat[,m$etamap$offsettheta,drop=FALSE])>1e-8,1,any)
    xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE]
    colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
    xmat <- xmat[!shouldoffset,,drop=FALSE]
    zy <- zy[!shouldoffset]
    wend <- wend[!shouldoffset]
    foffset <- foffset[!shouldoffset]
#    theta.offset <- theta.offset[!m$etamap$offsettheta]
  }else{
    foffset <- rep(0, length=length(zy))
    theta.offset <- rep(0, length=Clist$nparam)
    if(Clist$nedges>0){
      theta.offset[1] <- log(Clist$nedges/(Clist$ndyads-Clist$nedges))
    }else{
      theta.offset[1] <- log(1/(Clist$ndyads-1))
    }
    names(theta.offset) <- m$coef.names
  }
#
# Sample if necessary
#
  if(nrow(xmat) > maxMPLEsamplesize){
   rsample <- sample((1:nrow(xmat))[zy==1], size=min(maxMPLEsamplesize,sum(zy)),
                     replace=FALSE)
   rsample <- c(rsample, 
     sample((1:nrow(xmat))[zy==0], size=min(maxMPLEsamplesize,sum(!zy)),
                     replace=FALSE) )
   tau <- sum(zy*wend)/sum(wend)
   xmat.full <- xmat
   zy.full <- zy
   foffset.full <- foffset
   zy <- zy[rsample]
   wend <- wend[rsample]
   wend <- tau*zy*wend/sum(zy*wend) +
           (1-tau)*(1-zy)*wend/sum((1-zy)*wend)
   wend <- wend*nrow(xmat)/sum(wend)
   xmat <- xmat[rsample,,drop=FALSE]
   foffset <- foffset[rsample]
  }else{
   xmat.full <- NULL
   zy.full <- NULL
   foffset.full <- NULL
  }

  list(xmat=xmat, zy=zy, foffset=foffset, wend=wend, numobs=round(sum(wend)),
       xmat.full=xmat.full, zy.full=zy.full, foffset.full=foffset.full,
       theta.offset=theta.offset, maxMPLEsamplesize=maxMPLEsamplesize)
}
