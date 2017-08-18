#  File R/ergm.pl.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###############################################################################
# The <ergm.pl> function prepares many of the components needed by <ergm.mple>
# for the regression rountines that are used to find the MPLE estimated ergm;
# this is largely done through <MPLE_wrapper.C>.
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
#   control       : a list of MCMC related parameters; recognized variables
#                      include:
#         samplesize : the number of networks to sample, which will inform the size
#                      of the returned 'xmat'
#         Clist.miss : see 'Clist.miss' above; some of the code uses this Clist.miss,
#                      some uses the one above, does this seem right?
#   MHproposal       : an MHproposal object, as returned by <ergm.getMHproposal>
#   verbose          : whether this and the C routines should be verbose (T or F);
#                      default=FALSE
#
# --RETURNED--
#   a list containing
#     xmat     : the compressed and possibly sampled matrix of change
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
                    control, MHproposal, ignore.offset=FALSE,
                    verbose=FALSE) {
  bip <- Clist$bipartite
  n <- Clist$n

  maxNumDyadTypes <- min(control$MPLE.max.dyad.types,
                         ifelse(bip>0, bip*(n-bip), 
                                ifelse(Clist$dir, n*(n-1), n*(n-1)/2)))
                        
  # May have to think harder about what maxNumDyadTypes should be if we 
  # implement a hash-table approach to compression.
  # *** don't forget, pass in tails first now, not heads
  z <- .C("MPLE_wrapper",
          as.integer(Clist$tails), as.integer(Clist$heads),
          as.integer(Clist$nedges),
          as.integer(FALSE),
          as.integer(c(Clist.miss$nedges,Clist.miss$tails,Clist.miss$heads)),
          as.integer(n), 
          as.integer(Clist$dir),     as.integer(bip),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring), as.character(Clist$snamestring),
          as.double(Clist$inputs),
          y = integer(maxNumDyadTypes),
          x = double(maxNumDyadTypes*Clist$nstats),
          weightsvector = integer(maxNumDyadTypes),
          as.integer(.Machine$integer.max), # maxDyads
          as.integer(maxNumDyadTypes),
          PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  if (verbose) {
    message(paste("MPLE covariate matrix has", sum(uvals), "rows."))
  }
  zy <- z$y[uvals]
  wend <- as.numeric(z$weightsvector[uvals])
  xmat <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  rm(z,uvals)

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(wend)<Clist$ndyads && mean(zy)<1/2){
    if(verbose) message("A sparse network with too many unique dyads encountered. Using case-control MPLE.")
    # Strip out the rows associated with ties.
    wend <- wend[zy==0]
    xmat <- xmat[zy==0,,drop=FALSE]
    zy <- zy[zy==0]

    ## Run a whitelist PL over all of the edges in the network.
    # Generate a random permutation of edges, in case we run out of space here as well.
    missing <- paste(Clist$tails,Clist$heads,sep=".") %in% paste(Clist.miss$tails,Clist.miss$heads,sep=".")
    ordering <- sample.int(Clist$nedges-sum(missing))
    z <- .C("MPLE_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads),
            as.integer(Clist$nedges),
            as.integer(TRUE),
            as.integer(c(Clist$nedges-sum(missing),Clist$tails[!missing][ordering],Clist$heads[!missing][ordering])),
            as.integer(n), 
            as.integer(Clist$dir),     as.integer(bip),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring),
            as.double(Clist$inputs),
            y = integer(maxNumDyadTypes),
            x = double(maxNumDyadTypes*Clist$nstats),
            weightsvector = integer(maxNumDyadTypes),
            as.integer(.Machine$integer.max), # maxDyads
            as.integer(maxNumDyadTypes),
            PACKAGE="ergm")
    uvals <- z$weightsvector!=0
    zy.e <- z$y[uvals]
    wend.e <- as.numeric(z$weightsvector[uvals])
    xmat.e <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
    colnames(xmat.e) <- m$coef.names
    rm(z,uvals)

    # Divvy up the sampling weight of the ties:
    wend.e.total <- (Clist$nedges-sum(missing))
    wend.e <- wend.e / sum(wend.e) * wend.e.total

    # Divvy up the sampling weight of the nonties:
    wend <- wend / sum(wend) * (Clist$ndyads-wend.e.total)

    zy <- c(zy,zy.e)
    wend <- c(wend, wend.e)
    xmat <- rbind(xmat, xmat.e)

    rm(zy.e, wend.e, xmat.e)
  }
  

  #
  # Adjust for the offset
  #

  if(any(m$etamap$offsettheta) && !ignore.offset){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    # Compute the offset's effect.
    foffset <- .multiply.with.inf(xmat[,m$etamap$offsetmap,drop=FALSE], 
                                  cbind(ergm.eta(theta.offset,m$etamap)[m$etamap$offsetmap]))
    
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
