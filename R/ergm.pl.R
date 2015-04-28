#  File R/ergm.pl.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
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
#   conddeg          : an indicator of whether the MPLE should be conditional
#                      on degree; non-NULL values indicate yes, NULL no;
#                      default=NULL 
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
                    conddeg=NULL, control, MHproposal,
                    verbose=FALSE) {
  bip <- Clist$bipartite
  n <- Clist$n

  maxNumDyadTypes <- min(control$MPLE.max.dyad.types,
                         ifelse(bip>0, bip*(n-bip), 
                                ifelse(Clist$dir, n*(n-1), n*(n-1)/2)))
                        
  # May have to think harder about what maxNumDyadTypes should be if we 
  # implement a hash-table approach to compression.
  if(is.null(conddeg)){
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
    cat(paste("MPLE covariate matrix has", sum(uvals), "rows.\n"))
  }
  zy <- z$y[uvals]
  wend <- as.numeric(z$weightsvector[uvals])
  xmat <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  rm(z,uvals)

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(wend)<Clist$ndyads && mean(zy)<1/2){
    if(verbose) cat("A sparse network with too many unique dyads encountered. Using case-control MPLE.\n")
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
    data <- list(conddeg=conddeg,Clist=Clist,
         MHproposal=MHproposal, eta0=eta0,
         control=control,maxedges=maxedges,verbose=verbose)
    simfn <- function(i, data){
     # *** don't forget, pass in tails first now, not heads    
     z <- .C("MPLEconddeg_wrapper",
            as.integer(1), as.integer(data$conddeg$Clist$nedges),
            as.integer(data$conddeg$Clist$tails), as.integer(data$conddeg$Clist$heads),
            as.integer(data$conddeg$Clist$n),
            as.integer(data$conddeg$Clist$dir), as.integer(data$conddeg$Clist$bipartite),
            as.integer(data$conddeg$Clist$nterms),
            as.character(data$conddeg$Clist$fnamestring),
            as.character(data$conddeg$Clist$snamestring),
            as.character(data$MHproposal$name), as.character(data$MHproposal$pkgname),
            as.double(data$conddeg$Clist$inputs), as.double(data$eta0),
            as.integer(data$control$MPLE.samplesize),
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
            status = integer(1),
            PACKAGE="ergm")
    # save the results
    z <- list(s=z$s, newnwtails=z$newnwtails, newnwheads=z$newnwheads)
    
    nedges <- z$newnwtails[1]
    statsmatrix <- matrix(z$s, nrow=data$control$MPLE.samplesize+1,
                          ncol=data$conddeg$Clist$nstats,
                          byrow = TRUE)
    colnames(statsmatrix) <- data$conddeg$m$coef.names
    xb <- ergm.sufftoprob(statsmatrix[-c(1,data$control$MPLE.samplesize+1),],compress=TRUE)
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
  # =======================
  # Helper function
  # A is a matrix. V is a column vector that may contain Infs
  # computes A %*% V, counting 0*Inf as 0
  # May be slow if there are many rows. Use C here?
  multiply.with.inf <- function(A,V) {
    cbind(sapply(seq_len(nrow(A)), function(i) sum(V * A[i,], na.rm=T)))
  }

  if(any(m$etamap$offsettheta)){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    # Compute the offset's effect.
    foffset <- multiply.with.inf(xmat[,m$etamap$offsettheta,drop=FALSE], 
                                 cbind(theta.offset[m$etamap$offsettheta])) 
    
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
