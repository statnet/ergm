ergm.pl<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    maxMPLEsamplesize=1e+6,
                    maxNumDyadTypes=1e+6,
                    conddeg=NULL, MCMCparams, MHproposal,
                    verbose=FALSE, compressflag=TRUE) {
  bip <- Clist$bipartite
  n <- Clist$n
# Determine whether any edges are missing.  If so, this will reduce the
# number of observed (numobs) and also "turn off" the missing edges by
# setting offset to 1 in the corresponding rows.
  if(Clist.miss$nedges>0){
    temp <- matrix(0,ncol=n,nrow=n)
    base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
    base <- base[base[, 2] > base[, 1], ]
    if(Clist.miss$dir){
      base <- cbind(base[,c(2,1)],base)
      base <- matrix(t(base),ncol=2,byrow=TRUE)
    }
    ubase <- base[,1] + n*base[,2]
    offset <- !is.na(match(ubase, Clist.miss$heads+Clist.miss$tails*n))
    offset <- 1*offset
    numobs <- Clist$ndyads - sum(offset)
  }else{
    numobs <- Clist$ndyads
    offset <- rep(0,numobs)
  }

  maxNumDyadTypes <- min(maxNumDyadTypes,
                         ifelse(bip>0, bip*(n-bip), 
                                ifelse(Clist$dir, n*(n-1), n*(n-1)/2)))
  # May have to think harder about what maxNumDyadTypes should be if we 
  # implement a hash-table approach to compression.
  if(is.null(conddeg)){
  z <- .C("MPLE_wrapper",
          as.integer(Clist$heads),    as.integer(Clist$tails),
          as.integer(Clist$nedges),   as.integer(Clist$maxpossibleedges),
          as.integer(n), 
          as.integer(Clist$dir),     as.integer(bip),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring), as.character(Clist$snamestring),
          as.double(Clist$inputs),
          y = integer(maxNumDyadTypes),
          x = double(maxNumDyadTypes*Clist$nstats),
          weightsvector = integer(maxNumDyadTypes),
          as.double(offset), compressedOffset=double(maxNumDyadTypes),
          as.integer(maxNumDyadTypes),
          as.integer(maxMPLEsamplesize),
          as.integer(compressflag),
          PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  if (verbose) {
    if (compressflag) { cat("Compressed ") }
    cat(paste("MPLE covariate matrix has", sum(uvals), "rows.\n"))
  }
  zy <- z$y[uvals]
  wend <- z$weightsvector[uvals]
  xmat <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  dmiss <- z$compressedOffset[uvals]
  rm(z,uvals)
  }else{
    if (verbose) {
      cat("Using the MPLE conditional on degree.\n")
    }
    # Conditional on degree version
    eta0 <- ergm.eta(rep(0,length(conddeg$m$coef.names)), conddeg$m$etamap)
    
    stats <- matrix(0,ncol=conddeg$Clist$nstats,nrow=MCMCparams$samplesize+1)
    MCMCparams$stats <- stats
    maxedges <- max(5000, conddeg$Clist$nedges)
    flush.console()
    z <- .C("MPLEconddeg_wrapper",
            as.integer(conddeg$Clist$heads), as.integer(conddeg$Clist$tails),
            as.integer(conddeg$Clist$nedges), as.integer(conddeg$Clist$maxpossibleedges), 
            as.integer(conddeg$Clist$n),
            as.integer(conddeg$Clist$dir), as.integer(conddeg$Clist$bipartite),
            as.integer(conddeg$Clist$nterms),
            as.character(conddeg$Clist$fnamestring),
            as.character(conddeg$Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$package),
            as.double(conddeg$Clist$inputs), as.double(eta0),
            as.integer(MCMCparams$samplesize+1),
            s = as.double(t(MCMCparams$stats)),
            as.integer(0), 
            as.integer(1),
            newnwheads = integer(maxedges),
            newnwtails = integer(maxedges),
            as.integer(verbose), as.integer(MHproposal$bd$attribs),
            as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
            as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
            as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
            as.integer(maxedges),
            as.integer(MCMCparams$Clist.miss$heads), as.integer(MCMCparams$Clist.miss$tails),
            as.integer(MCMCparams$Clist.miss$nedges),
            PACKAGE="ergm")
    # save the results
    z <- list(s=z$s, newnwheads=z$newnwheads, newnwtails=z$newnwtails)
    
    nedges <- z$newnwheads[1]
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize+1,
                          ncol=conddeg$Clist$nstats,
                          byrow = TRUE)
    colnames(statsmatrix) <- conddeg$m$coef.names
    # xb <- apply(statsmatrix,2,diff)
    xb <- statsmatrix[-1,]
    zy <- round(xb[,1]-1)
    xmat <- xb[,-1,drop=FALSE]
    xmat[zy==1,] <- -xmat[zy==1,]
    wend <- zy-zy+1
    #
    #foffset <- NULL
    #foffset.full <- NULL
    #xmat.full <- xmat
    #zy.full <- zy
  }

  #
  # Adjust for the offset
  #
  if(any(m$etamap$offsettheta)){
    if(is.null(theta.offset)){
      theta.offset <- rep(0, length=Clist$nstats)
      names(theta.offset) <- m$coef.names
      theta.offset[m$etamap$offsettheta] <- -Inf
#     theta.offset[m$etamap$offsettheta] <- -10000
    }
# Commenting out recent version of this section that does not work;
#   foffset <- xmat[,m$etamap$offsettheta,drop=FALSE]%*%theta.offset[m$etamap$offsettheta]
#  Returning to the version from CRAN ergm v. 2.1:

    foffset <- xmat[,!m$etamap$offsettheta,drop=FALSE]%*%theta.offset[!m$etamap$offsettheta]
    shouldoffset <- apply(abs(xmat[,m$etamap$offsettheta,drop=FALSE])>1e-8,1,any)
    xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE]
    colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
    xmat <- xmat[!shouldoffset,,drop=FALSE]
    zy <- zy[!shouldoffset]
    wend <- wend[!shouldoffset]
    foffset <- foffset[!shouldoffset]
#   theta.offset <- theta.offset[!m$etamap$offsettheta]
  }else{
    foffset <- rep(0, length=length(zy))
    theta.offset <- rep(0, length=Clist$nstats)
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
    if(is.null(colnames(xmat))){colnames(xmat) <- m$coef.names}
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
