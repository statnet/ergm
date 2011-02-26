###############################################################################
# The <ergm.pl> function prepares many of the components needed by <ergm.mple>
# for the regression rountines that are used to find the MPLE estimated ergm;
# this is largely done through <MPLE_wrapper.C> or <MPLEconddeg_wrapper>
#
# --PARAMETERS--
#   Clist            : a list of parameters used for fitting and returned
#                      by <ergm.Cprepare>
#   Clist.miss       : the corresponding 'Clist' for the network of missing
#                      edges returned by <ergm.design>
#   m                : the model, as returned by <ergm.getmodel>
#   theta.offset     : a logical vector specifying which of the model
#                      coefficients are offset, i.e. fixed
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the
#                      pseudolikelihood computation; default=1e6
#   maxNumDyadTypes  : the maximum number of unique pseudolikelihood
#                      change statistics to be allowed if 'compressflag'=TRUE;
#                      if this is less than the networks maximum count of
#                      unique change stats, these will be sampled to attain
#                      the correct size; default=1e6
#   conddeg          : an indicator of whether the MPLE should be conditional
#                      on degree; non-NULL values indicate yes, NULL no;
#                      default=NULL 
#   MCMCparams       : a list of MCMC related parameters; recognized variables
#                      include:
#         samplesize : the number of networks to sample, which will inform the size
#                      of the returned 'xmat'
#         Clist.miss : see 'Clist.miss' above; some of the code uses this Clist.miss,
#                      some uses the one above, does this seem right?
#   MHproposal       : an MHproposal object, as returned by <ergm.getMHproposal>
#   verbose          : whether this and the C routines should be verbose (T or F);
#                      default=FALSE
#   compressflag     : whether to compress the design matrix of change stats by
#                      tabulating the unique rows (T or F); default=TRUE
#
#
# --RETURNED--
#   a list containing
#     xmat     : the possibly compressed and possibly sampled matrix of change
#                statistics
#     zy       : the corresponding vector of responses, i.e. tie values
#     foffset  : ??
#     wend     : the vector of weights for 'xmat' and 'zy'
#     numobs   : the number of dyads 
#     xmat.full: the 'xmat' before sampling; if no sampling is needed, this
#                is NULL
#     zy.full  : the 'zy' before  sampling; if no sampling is needed, this
#                is NULL
#     foffset.full     : ??
#     theta.offset     : a numeric vector whose ith entry tells whether the
#                        the ith curved coefficient?? was offset/fixed; -Inf
#                        implies the coefficient was fixed, 0 otherwise; if
#                        the model hasn't any curved terms, the first entry
#                        of this vector is one of
#                           log(Clist$nedges/(Clist$ndyads-Clist$nedges))
#                           log(1/(Clist$ndyads-1))
#                        depending on 'Clist$nedges'
#     maxMPLEsamplesize: the 'maxMPLEsamplesize' inputted to <ergm.pl>
#    
###############################################################################

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
    # *** don't forget, pass in tails first now, not heads    
    z <- .C("MPLEconddeg_wrapper",
            as.integer(conddeg$Clist$tails), as.integer(conddeg$Clist$heads),
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
            newnwtails = integer(maxedges),
            newnwheads = integer(maxedges),
            as.integer(verbose), as.integer(MHproposal$bd$attribs),
            as.integer(MHproposal$bd$maxout), as.integer(MHproposal$bd$maxin),
            as.integer(MHproposal$bd$minout), as.integer(MHproposal$bd$minin),
            as.integer(MHproposal$bd$condAllDegExact), as.integer(length(MHproposal$bd$attribs)),
            as.integer(maxedges),
            as.integer(MCMCparams$Clist.miss$tails), as.integer(MCMCparams$Clist.miss$heads),
            as.integer(MCMCparams$Clist.miss$nedges),
            PACKAGE="ergm")
    # save the results
    z <- list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads)
    
    nedges <- z$newnwtails[1]
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
    #foffset <- xmat[,m$etamap$offsettheta,drop=FALSE]%*%theta.offset[m$etamap$offsettheta]
    #foffset[is.nan(foffset)] <- 0 # zero times +-Inf should be zero in this context
    #foffset <- xmat[,m$etamap$offsettheta,drop=FALSE]%*%theta.offset[m$etamap$offsettheta]
    #xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE]
    #colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
    
    # Returning to the version from CRAN ergm v. 2.1:
    foffset <- xmat[,!m$etamap$offsettheta,drop=FALSE] %*%
               theta.offset[!m$etamap$offsettheta]
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
