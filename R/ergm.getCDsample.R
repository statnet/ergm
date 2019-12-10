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
                             response=NULL, verbose=FALSE,..., eta=ergm.eta(theta, model$etamap), stats0=NULL){
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.network(nw) || is.ergm_state(nw)) nw <- list(nw)
  nws <- rep(nw, length.out=nthreads(control))
  stats <- if(is.network(nws[[1]])){
    NVL(stats0) <- numeric(length(eta))
    if(is.numeric(stats0)) stats0 <- list(stats0)
    stats0 <- rep(stats0, length.out=length(nws))
    
    states <- mapply(ergm_state, nws, stats0, MoreArgs=list(model=model, response=response), SIMPLIFY=FALSE)
  }else nws

  control.parallel <- control
  control.parallel$CD.samplesize <- NVL3(control$CD.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  doruns <- function(samplesize=NULL){
    if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control), ergm_CD_slave,
                                  state=states, MoreArgs=list(model=model,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...,samplesize=samplesize))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_CD_slave(state=states[[1]], samplesize=samplesize,model=model,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  outl <- doruns()
  for(i in seq_along(outl)){
    outl[[i]]$s <- mcmc(outl[[i]]$s)
  }
  
  if(control.parallel$MCMC.runtime.traceplot){
      esteq <- lapply.mcmc.list(lapply(outl, function(out)
                      NVL3(theta, ergm.estfun(out$s, ., model), out$s[,!model$offsetmap,drop=FALSE])
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

ergm_CD_slave <- function(state, model, proposal,eta,control,verbose,..., samplesize=NULL, stats0=numeric(nparam(model, canonical=TRUE))){
  stats <- stats0
  
  if(is.null(samplesize)) samplesize <- control$CD.samplesize
  
  z <- if(!is.valued(state))
         .Call("CD_wrapper",
               state, model, proposal,
               # MCMC settings
               as.double(.deinf(eta)),
               as.integer(samplesize),
               as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
               as.integer(verbose),
               PACKAGE="ergm")
       else
         .Call("WtCD_wrapper",
               state, model, proposal,
               # MCMC settings
               as.double(.deinf(eta)),
               as.integer(samplesize),
               as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
               as.integer(verbose),
               PACKAGE="ergm")

  z$s <- matrix(z$s+stats, ncol=nparam(model,canonical=TRUE), byrow = TRUE)

  z
}
