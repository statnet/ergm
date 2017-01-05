#  File R/ergm.getMCMCsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#==============================================================================
# This file contains the 2 following functions for getting an MCMC sample
#      <ergm.getMCMCsample>
#      <ergm.mcmcslave>
#==============================================================================




#########################################################################################
# The <ergm.getMCMCsample> function samples networks using an MCMC algorithm via
# <MCMC_wrapper.C>. Unlike its <ergm.getMCMCsample> counterpart, this function is
# caple of running in multiple threads.  Note that the returned stats will be relative to
# the original network, i.e., the calling function must shift the statistics if required. 
# The calling function must also attach column names to the statistics matrix if required.
#
# --PARAMETERS--
#   nw        :  a network object
#   model     :  a model for the given 'nw' as returned by <ergm.getmodel>
#   MHproposal:  a list of the parameters needed for Metropolis-Hastings proposals and
#                the result of calling <MHproposal>
#   eta0      :  the initial eta coefficients
#   verbose   :  whether the C functions should be verbose; default=FALSE
#   control:  list of MCMC tuning parameters; those recognized include
#       parallel    : the number of threads in which to run the sampling
#       packagenames: names of packages; this is only relevant if "ergm" is given
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, control, verbose should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the sample as a list containing:
#     statsmatrix:  the stats matrix for the sampled networks, RELATIVE TO THE ORIGINAL
#                   NETWORK!
#     newnetwork :  the edgelist of the final sampled network
#     nedges     :  the number of edges in the 'newnetwork' ??
#
#########################################################################################

ergm.getMCMCsample <- function(nw, model, MHproposal, eta0, control, 
                                        verbose, response=NULL, ...) {
  nthreads <- max(
    if(inherits(control$parallel,"cluster")) nrow(summary(control$parallel))
    else control$parallel,
    1)

  cl <- if(!is.numeric(control$parallel) || control$parallel!=0){
    ergm.getCluster(control, verbose)
  }else NULL
  
  if(is.network(nw)) nw <- list(nw)
  nws <- rep(nw, length.out=nthreads)
  
  Clists <- lapply(nws, ergm::ergm.Cprepare, model, response=response)

  control.parallel <- control
  if(!is.null(control$MCMC.samplesize)) control.parallel$MCMC.samplesize <- ceiling(control$MCMC.samplesize / nthreads)

  flush.console()

  doruns <- function(prev.runs=rep(list(NULL),nthreads), burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL){
    if(!is.null(cl)) clusterMap(cl,ergm.mcmcslave,
                                  Clist=Clists, prev.run=prev.runs, MoreArgs=list(MHproposal=MHproposal,eta0=eta0,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval,maxedges=maxedges))
    else list(ergm.mcmcslave(Clist=Clists[[1]], prev.run=prev.runs[[1]],burnin=burnin,samplesize=samplesize,interval=interval,maxedges=maxedges,MHproposal=MHproposal,eta0=eta0,control=control.parallel,verbose=verbose,...))
  }
  
  if(!is.null(control.parallel$MCMC.effectiveSize)){
    if(verbose) cat("Beginning adaptive MCMC...\n")

    howmuchmore <- function(target.ess, current.ss, current.ess, current.burnin){
      (target.ess-current.ess)*(current.ss-current.burnin)/current.ess
    }
    
    interval <- control.parallel$MCMC.interval
    meS <- list(burnin=0,eS=control.parallel$MCMC.effectiveSize)
    outl <- rep(list(NULL),nthreads)
    for(mcrun in seq_len(control.parallel$MCMC.effectiveSize.maxruns)){
      if(mcrun==1){
        samplesize <- control.parallel$MCMC.samplesize
        if(verbose)
          cat("First run: running each chain forward by",samplesize, "steps with interval", interval, ".\n")
      }else{
        if(meS$eS<1){
          samplesize <- control.parallel$MCMC.samplesize
          if(verbose)
            cat("Insufficient ESS to determine the number of steps remaining: running forward by",samplesize, "steps with interval", interval, ".\n")
        }else{
          pred.ss <- howmuchmore(control.parallel$MCMC.effectiveSize, NVL(nrow(outl[[1]]$s),0), meS$eS, meS$burnin)
          damp.ss <- pred.ss*(meS$eS/(control.parallel$MCMC.effectiveSize.damp+meS$eS))+control.parallel$MCMC.samplesize*(1-meS$eS/(control.parallel$MCMC.effectiveSize.damp+meS$eS))
          samplesize <- round(damp.ss)
          if(verbose) cat("Predicted additional sample size:",pred.ss, "dampened to",damp.ss, ", so running", samplesize, "steps forward.\n")
        }
      }
        
      outl<-doruns(prev.runs=outl,
                  burnin = interval, # I.e., skip that much before the first draw.
                  samplesize = samplesize,
                  interval = interval,
                  maxedges = if(!is.null(outl[[1]])) max(sapply(outl,"[[","maxedges")) else NULL
                  )
      
      # Stop if something went wrong.
      if(any(sapply(outl,"[[","status")!=0)) break
      
      while(nrow(outl[[1]]$s)-meS$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)) outl[[i]]$s <- outl[[i]]$s[seq_len(floor(nrow(outl[[i]]$s)/2))*2,,drop=FALSE]
        interval <- interval*2
        if(verbose) cat("Increasing thinning to",interval,".\n")
      }
      
      esteq <- lapply(outl, function(out)
                      if(all(c("theta","etamap") %in% names(list(...)))) .ergm.esteq(list(...)$theta, list(etamap=list(...)$etamap), out$s)
                      else out$s[,Clists[[1]]$diagnosable,drop=FALSE]
                      )
      
      meS <- .max.effectiveSize(esteq, npts=control$MCMC.effectiveSize.points, base=control$MCMC.effectiveSize.base, ar.order=control$MCMC.effectiveSize.order)
      if(verbose) cat("Maximum Harmonic mean ESS of",meS$eS,"attained with burn-in of", round(meS$b/nrow(outl[[1]]$s)*100,2),"%.\n")

      if(control.parallel$MCMC.runtime.traceplot){
        for (i in seq_along(esteq)) colnames(esteq[[i]]) <- names(list(...)$theta)
        plot(coda::as.mcmc.list(lapply(lapply(esteq, coda::mcmc), window, thin=max(1,floor(nrow(esteq)/1000))))
             ,ask=FALSE,smooth=TRUE,density=FALSE)
      }

      if(meS$eS>=control.parallel$MCMC.effectiveSize){
        if(verbose) cat("Target ESS achieved. Returning.\n")      
        break
      }
    }
    
    if(meS$eS<control.parallel$MCMC.effectiveSize)
      warning("Unable to reach target effective size in iterations alotted.")

    for(i in seq_along(outl)){
      if(meS$burnin) outl[[i]]$s <- outl[[i]]$s[-seq_len(meS$burnin),,drop=FALSE]
      outl[[i]]$s <- coda::mcmc(outl[[i]]$s, (meS$burnin+1)*interval, thin=interval)
      outl[[i]]$final.interval <- interval
    }
  }else{
    outl <- doruns()
    for(i in seq_along(outl)){
      outl[[i]]$s <- coda::mcmc(outl[[i]]$s, control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
    }
    
    if(control.parallel$MCMC.runtime.traceplot){
      esteq <- lapply(outl, function(out)
        if(all(c("theta","etamap") %in% names(list(...)))) .ergm.esteq(list(...)$theta, list(etamap=list(...)$etamap), out$s)
        else out$s[,Clists[[1]]$diagnosable,drop=FALSE]
      )
      for (i in seq_along(esteq)) colnames(esteq[[i]]) <- names(list(...)$theta)
      plot(coda::as.mcmc.list(lapply(lapply(esteq, coda::mcmc), window, thin=max(1,floor(nrow(esteq)/1000))))
           ,ask=FALSE,smooth=TRUE,density=FALSE)
    }
  }

  #
  #   Process the results
  #
  statsmatrices <- list()
  newnetworks <- list()
  final.interval <- c()
  for(i in (1:nthreads)){
    z <- outl[[i]]
    
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control.parallel$MCMC.max.maxedges
      return(list(status=z$status))
    }
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }
    
    statsmatrices[[i]] <- z$s
    
    newnetworks[[i]]<-newnw.extract(nws[[i]],z,
                                    response=response,output=control.parallel$network.output)
    final.interval <- c(final.interval, z$final.interval)
  }
  newnetwork<-newnetworks[[1]]
  
  
  ergm.stopCluster(cl)

  statsmatrix <- do.call(rbind,statsmatrices)
  colnames(statsmatrix) <- model$coef.names

  if(verbose){cat("Sample size =",nrow(statsmatrix),"by",
                  control.parallel$MCMC.samplesize,"\n")}
  
  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, statsmatrices=statsmatrices, newnetwork=newnetwork, newnetworks=newnetworks, status=0, final.interval=final.interval)

}



###############################################################################
# The <ergm.mcmcslave> function is that which the slaves will call to perform
# a validation on the mcmc equal to their slave number. It also returns an
# MCMC sample.
#
# --PARAMETERS--
#   Clist     : the list of parameters returned by <ergm.Cprepare>
#   MHproposal: the MHproposal list as returned by <getMHproposal>
#   eta0      : the canonical eta parameters
#   control: a list of parameters for controlling the MCMC algorithm;
#               recognized components include:
#       samplesize  : the number of networks to be sampled
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#   verbose   : whether the C code should be verbose (T or F)
#   ...       : optional arguments
#
# --RETURNED--
#   the MCMC sample as a list of the following:
#     s         : the statsmatrix
#     newnwtails: the vector of tails for the new network- is this the final
#                 network sampled? - is this the original nw if 'maxedges' is 0
#     newnwheads: the vector of heads for the new network - same q's
#
###############################################################################

ergm.mcmcslave <- function(Clist,MHproposal,eta0,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL) {

  numnetworks <- 0

  if(is.null(prev.run)){ # Start from Clist
    nedges <- c(Clist$nedges,0,0)
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
    stats <- rep(0, Clist$nstats)
  }else{ # Pick up where we left off
    nedges <- prev.run$newnwtails[1]
    tails <- prev.run$newnwtails[2:(nedges+1)]
    heads <- prev.run$newnwheads[2:(nedges+1)]
    weights <- prev.run$newnwweights[2:(nedges+1)]
    nedges <- c(nedges,0,0)
    stats <- prev.run$s[nrow(prev.run$s),]
  }
  
  if(is.null(burnin)) burnin <- control$MCMC.burnin
  if(is.null(samplesize)) samplesize <- control$MCMC.samplesize
  if(is.null(interval)) interval <- control$MCMC.interval
  if(is.null(maxedges)) maxedges <- control$MCMC.init.maxedges
  
  repeat{
    if(is.null(Clist$weights)){
      z <- .C("MCMC_wrapper",
              as.integer(numnetworks), as.integer(nedges),
              as.integer(tails), as.integer(heads),
              as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(MHproposal$name), as.character(MHproposal$pkgname),
              as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
              as.integer(samplesize),
              s = as.double(rep(stats, samplesize)),
              as.integer(burnin), 
              as.integer(interval),
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
              as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
              as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
              as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
              as.integer(maxedges),
              status = integer(1),
              PACKAGE="ergm")
      
        # save the results (note that if prev.run is NULL, c(NULL$s,z$s)==z$s.
      z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE),
                newnwtails=z$newnwtails, newnwheads=z$newnwheads, status=z$status, maxedges=maxedges)
    }else{
      z <- .C("WtMCMC_wrapper",
              as.integer(length(nedges)), as.integer(nedges),
              as.integer(tails), as.integer(heads), as.double(weights),
              as.integer(Clist$n),
              as.integer(Clist$dir), as.integer(Clist$bipartite),
              as.integer(Clist$nterms),
              as.character(Clist$fnamestring),
              as.character(Clist$snamestring),
              as.character(MHproposal$name), as.character(MHproposal$pkgname),
              as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
              as.integer(samplesize),
              s = as.double(rep(stats, samplesize)),
              as.integer(burnin), 
              as.integer(interval),
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              newnwweights = double(maxedges),
              as.integer(verbose), 
              as.integer(maxedges),
              status = integer(1),
              PACKAGE="ergm")
      # save the results
      z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE),
              newnwtails=z$newnwtails, newnwheads=z$newnwheads, newnwweights=z$newnwweights, status=z$status, maxedges=maxedges)
    }
    
    z$s <- rbind(prev.run$s,z$s)
    colnames(z$s) <- names(Clist$diagnosable)
    
    if(z$status!=1) return(z) # Handle everything except for MCMC_TOO_MANY_EDGES elsewhere.
    
    # The following is only executed (and the loop continued) if too many edges.
    maxedges <- maxedges * 10
    if(!is.null(control$MCMC.max.maxedges)){
      if(maxedges == control$MCMC.max.maxedges*10) # True iff the previous maxedges exactly equaled control$MCMC.max.maxedges and that was too small.
        return(z) # This will kick the too many edges problem upstairs, so to speak.
      maxedges <- min(maxedges, control$MCMC.max.maxedges)
    }
  }
}


.max.effectiveSize <- function(x, npts, base, ar.order=0){
  if(!is.list(x)) x <- list(x)
  es <- function(b){
    if(b>0) x <- lapply(lapply(x, "[", -seq_len(b),,drop=FALSE),coda::mcmc)
    effSizes <- if(ar.order) .fast.effectiveSize(as.matrix(coda::as.mcmc.list(x), ar.order=ar.order))
                else effectiveSize(as.matrix(coda::as.mcmc.list(x)))
    mean.fn <- function(x) x^(-1)
    mean.ifn <- function(x) x^(-1)
    mean.ifn(mean(mean.fn(effSizes)))
  }

  pts <- sort(round(base^seq_len(npts)*nrow(x[[1]])))
  ess <- sapply(pts, es)

  list(burnin=pts[which.max(ess)], eS=max(ess))
}

.fast.effectiveSize <- function(x, ar.order=1){
  if (coda::is.mcmc.list(x)){
    ess <- do.call(rbind, lapply(x, .fast.effectiveSize, ar.order=ar.order))
    ans <- apply(ess, 2, base::sum)
  } else {
    x <- coda::as.mcmc(x)
    x <- as.matrix(x)
    spec <- .fast.spectrum0.ar(x, ar.order=ar.order)$spec
    ans <- ifelse(spec == 0, 0, nrow(x) * apply(x, 2, stats::var)/spec)
    }
    return(ans)
}
.fast.spectrum0.ar <- function (x, ar.order=1){
    x <- as.matrix(x)
    v0 <- order <- numeric(ncol(x))
    names(v0) <- names(order) <- colnames(x)
    z <- 1:nrow(x)
    for (i in 1:ncol(x)) {
      ar.out <- ar(x[, i], aic = FALSE, order.max=ar.order)
      v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
      order[i] <- ar.out$order
    }
    return(list(spec = v0, order = order))
}
