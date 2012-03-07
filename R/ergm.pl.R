#  File ergm/R/ergm.pl.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <ergm.pl> function prepares many of the components needed by <ergm.mple>
# for the regression rountines that are used to find the MPLE estimated ergm;
# this is largely done through <MPLE_wrapper.C> or <MPLEconddeg_wrapper>
###############################################################################

ergm.pl<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    maxMPLEsamplesize=1e+6,
                    maxNumDyadTypes=1e+6,
                    conddeg=NULL, control, MHproposal,
                    verbose=FALSE, compressflag=TRUE) {
  bip <- Clist$bipartite
  n <- Clist$n

  # Determine whether any edges are missing.  If so, this will reduce the
  # number of observed (numobs) and also "turn off" the missing edges by
  # setting offset to 1 in the corresponding rows.
  if(Clist.miss$nedges>0){
    nrows <- ifelse(bip>0, bip, n)
    ncolumns <- ifelse(bip>0, n-bip, n)
    base <- cbind(rep(1:nrows, rep(ncolumns, nrows)),
                  rep((1+bip):(ncolumns+bip), nrows))
    base <- base[base[, 2] > base[, 1], ]
    # At this point, the rows of base are in dictionary order
    if(Clist.miss$dir){
      # If we get here, then interleave the transposed rows, as in
      # (1,2), (2,1),   (1,3), (3,1),   (1,4), (4,1),  etc.
      # This is necessary because it's the order in which the
      # MPLE wrapper expects them to be entered, which is only
      # important when there are missing data because the "offset" 
      # vector in that case must index the correct rows.
      # (The C code reconstructs the base matrix from scratch.)
      base <- cbind(base, base[,c(2,1)])
      base <- matrix(t(base),ncol=2,byrow=TRUE)
    }
    ubase <- base[,1] + n*base[,2]
    offset <- !is.na(match(ubase, Clist.miss$tails+Clist.miss$heads*n))
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
  # *** don't forget, pass in tails first now, not heads
  z <- .C("MPLE_wrapper",
          as.integer(Clist$tails),    as.integer(Clist$heads),
          as.integer(Clist$nedges),   
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
    
    maxedges <- max(5000, conddeg$Clist$nedges)
    nsim <- 1
    if(control$MPLE.samplesize > 50000){
     nsim <- ceiling(control$MPLE.samplesize / 50000)
     control$MPLE.samplesize <- 50000

     cl <- ergm.getCluster(control, verbose)
    }
    flush.console()
#
    control$stats <- matrix(0,ncol=conddeg$Clist$nstats,nrow=control$MPLE.samplesize+1)
    data <- list(conddeg=conddeg,Clist=Clist, MHproposal=MHproposal, eta0=eta0,
         control=control,maxedges=maxedges,verbose=verbose)
    simfn <- function(i, data){
     # *** don't forget, pass in tails first now, not heads    
     z <- .C("MPLEconddeg_wrapper",
            as.integer(data$conddeg$Clist$tails), as.integer(data$conddeg$Clist$heads),
            as.integer(data$conddeg$Clist$nedges), as.integer(data$conddeg$Clist$maxpossibleedges), 
            as.integer(data$conddeg$Clist$n),
            as.integer(data$conddeg$Clist$dir), as.integer(data$conddeg$Clist$bipartite),
            as.integer(data$conddeg$Clist$nterms),
            as.character(data$conddeg$Clist$fnamestring),
            as.character(data$conddeg$Clist$snamestring),
            as.character(data$MHproposal$name), as.character(data$MHproposal$package),
            as.double(data$conddeg$Clist$inputs), as.double(data$eta0),
            as.integer(data$control$MPLE.samplesize+1),
            s = as.double(t(data$control$stats)),
            as.integer(0), 
            as.integer(1),
            newnwtails = integer(data$maxedges),
            newnwheads = integer(data$maxedges),
            as.integer(data$verbose), as.integer(data$MHproposal$arguments$constraints$bd$attribs),
            as.integer(data$MHproposal$arguments$constraints$bd$maxout), as.integer(data$MHproposal$arguments$constraints$bd$maxin),
            as.integer(data$MHproposal$arguments$constraints$bd$minout), as.integer(data$MHproposal$arguments$constraints$bd$minin),
            as.integer(data$MHproposal$arguments$constraints$bd$condAllDegExact),
as.integer(length(data$MHproposal$arguments$constraints$bd$attribs)),
            as.integer(data$maxedges),
            as.integer(data$control$Clist.miss$tails), as.integer(data$control$Clist.miss$heads),
            as.integer(data$control$Clist.miss$nedges),
            PACKAGE="ergm")
    # save the results
    z <- list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads)
    
    nedges <- z$newnwtails[1]
    statsmatrix <- matrix(z$s, nrow=data$control$MPLE.samplesize+1,
                          ncol=data$conddeg$Clist$nstats,
                          byrow = TRUE)
    colnames(statsmatrix) <- data$conddeg$m$coef.names
    xb <- ergm.sufftoprob(statsmatrix[-1,],compress=TRUE)
# if (verbose) {cat("Finished compression.\n")}
    xmat <- xb[,-c(1,ncol(xb)),drop=FALSE]
    wend <- xb[,ncol(xb)]
#   xb <- statsmatrix[-1,]
    zy <- round(xb[,1]-1)
    xmat[zy==1,] <- -xmat[zy==1,]
#   wend <- zy-zy+1
    return(list(zy=zy,xmat=xmat,wend=wend))
    }
#
    if(nsim >1){
      outlist <- clusterApplyLB(cl, as.list(1:nsim), simfn, data)
#
#     Process the results
#
      zy <- NULL
      xmat <- NULL
      wend <- NULL
      for(i in (1:nsim)){
       z <- outlist[[i]]
       zy <- c(zy, z$zy)
       xmat <- rbind(xmat, z$xmat)
       wend <- c(wend, z$wend)
      }
      rm(outlist)
      ergm.stopCluster(cl)
     }else{
      z <- simfn(1, data)
      zy <- z$zy
      xmat <- z$xmat
      wend <- z$wend
      rm(z)
     }
  }

  #
  # Adjust for the offset
  #
  if(any(m$etamap$offsettheta)){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    foffset <- xmat[,m$etamap$offsettheta,drop=FALSE] %*% cbind(theta.offset[m$etamap$offsettheta]) # Compute the offset's effect.
    foffset[is.nan(foffset)] <- 0 # 0*Inf==0 in this case.
    # Remove offset covariate columns.
    xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE] 
    colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
    # Now, iff a row's offset effect is infinite, then it carries no
    # further information whatsoever, so it should be dropped.
    xmat <- xmat[is.finite(foffset),,drop=FALSE]
    zy <- zy[is.finite(foffset)]
    wend <- wend[is.finite(foffset)]
    foffset <- foffset[is.finite(foffset)]
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
