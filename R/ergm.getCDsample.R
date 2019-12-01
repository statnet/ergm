#  File R/ergm.getCDsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

ergm_CD_sample <- function(nw, model, proposal, control, theta=NULL, 
                             response=NULL, verbose=FALSE,..., eta=ergm.eta(theta, model$etamap), stats0=NULL) {
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.network(nw) || is.pending_update_network(nw)) nw <- list(nw)
  nws <- rep(nw, length.out=nthreads(control))
  NVL(stats0) <- numeric(length(eta))
  if(is.numeric(stats0)) stats0 <- list(stats0)
  stats0 <- rep(stats0, length.out=length(nws))

  Clists <- lapply(nws, ergm::ergm.Cprepare, model, response=response)

  control.parallel <- control
  control.parallel$CD.samplesize <- NVL3(control$CD.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  doruns <- function(samplesize=NULL){
    if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control), ergm_CD_slave,
                                  Clist=Clists, stats0=stats0, MoreArgs=list(proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...,samplesize=samplesize))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_CD_slave(Clist=Clists[[1]], stats0=stats0[[1]], samplesize=samplesize,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  outl <- doruns()
  for(i in seq_along(outl)){
    outl[[i]]$s <- mcmc(outl[[i]]$s)
  }
  
  if(control.parallel$MCMC.runtime.traceplot){
      esteq <- lapply.mcmc.list(lapply(outl, function(out)
                      NVL3(theta, ergm.estfun(out$s, ., model), out$s[,Clists[[1]]$diagnosable,drop=FALSE])
                      ), mcmc, start=1)
        plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
             ,ask=FALSE,smooth=TRUE,density=FALSE)
  }

  #
  #   Process the results
  #
  statsmatrices <- list()
  newnetworks <- list()
  for(i in (1:nthreads(control))){
    z <- outl[[i]]
        
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
    }
    
    statsmatrices[[i]] <- z$s
  }
  
  stats <- as.mcmc.list(statsmatrices)
  if(verbose){message("Sample size = ",niter(stats)*nchain(stats)," by ",
                  niter(stats),".")}
  
  list(stats = stats, networks=newnetworks, status=0)
}

ergm_CD_slave <- function(Clist,proposal,eta,control,verbose,..., samplesize=NULL, stats0=numeric(Clist$nstats)) {
    nedges <- Clist$nedges
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
    stats <- stats0
  
  if(is.null(samplesize)) samplesize <- control$CD.samplesize
  
  z <-
    if(is.null(Clist$weights)){
      .Call("CD_wrapper",
            # Network settings
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            # Model settings
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            # Proposal settings
            as.character(proposal$name), as.character(proposal$pkgname),
            as.integer(proposal$arguments$constraints$bd$attribs),
            as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
            as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
            as.integer(proposal$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal$arguments$constraints$bd$attribs)),
            # Numeric vector inputs
            as.double(c(Clist$inputs,Clist$slots.extra.aux,proposal$inputs)),
            # Network state
            as.integer(nedges),
            as.integer(tails), as.integer(heads),
            # MCMC settings
            as.double(.deinf(eta)),
            as.integer(samplesize),
            as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
            as.integer(verbose),
            PACKAGE="ergm")
    }else{
      .Call("WtCD_wrapper",
            # Network settings
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            # Model settings
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            # Proposal settings
            as.character(proposal$name), as.character(proposal$pkgname),
            # Numeric inputs
            as.double(c(Clist$inputs,Clist$slots.extra.aux,proposal$inputs)),
            # Network state
            as.integer(nedges),
            as.integer(tails), as.integer(heads), as.double(weights),
            # MCMC settings
            as.double(.deinf(eta)),
            as.integer(samplesize),
            as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
            as.integer(verbose), 
            PACKAGE="ergm")
    }
    # save the results
  z<-list(s=matrix(z[[2]]+stats, ncol=Clist$nstats, byrow = TRUE),
          status=z[[1]])
}
