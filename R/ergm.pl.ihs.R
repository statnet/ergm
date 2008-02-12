ergm.pl.ihs<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    maxMPLEsamplesize=100000,
                    maxNumDyadTypes=100000,
                    verbose=FALSE) {
  if(Clist.miss$nedges>0){
    temp <- matrix(0,ncol=Clist$n,nrow=Clist$n)
    base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
    base <- base[base[, 2] > base[, 1], ]
    if(Clist.miss$dir){
      base <- cbind(base[,c(2,1)],base)
      base <- matrix(t(base),ncol=2,byrow=TRUE)
    }
    ubase <- base[,1] + Clist$n*base[,2]
    offset <- !is.na(match(ubase, Clist.miss$heads+Clist.miss$tails*Clist$n))
    offset <- 1*offset
    numobs <- Clist$ndyads - sum(offset)
  }else{
    offset <- rep(0,Clist$ndyads)
    numobs <- Clist$ndyads
  }
  z <- .C("MPLE_wrapper",
          as.integer(Clist$heads),    as.integer(Clist$tails),
          as.integer(Clist$nedges),   as.integer(Clist$maxpossibleedges),
          as.integer(Clist$n), 
          as.integer(Clist$dir),     as.integer(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring), as.character(Clist$snamestring),
          as.double(Clist$inputs),
          y = integer(maxNumDyadTypes),
          x = double(maxNumDyadTypes*Clist$nparam),
          weightsvector = integer(maxNumDyadTypes),
          as.double(offset), compressedOffset=double(maxNumDyadTypes),
          as.integer(maxNumDyadTypes),
          as.integer(maxMPLEsamplesize),
          PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  zy <- z$y[uvals]
  wend <- z$weightsvector[uvals]
  xmat <- matrix(z$x, ncol=Clist$nparam, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  dmiss <- z$compressedOffset[uvals]
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
  
  if(Clist.miss$nedges>0){
    xmat <- xmat[dmiss==0,,drop=FALSE]
    zy <- zy[dmiss==0]
    wend <- wend[dmiss==0]
    foffset <- foffset[dmiss==0]
    colnames(xmat) <- m$coef.names
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
