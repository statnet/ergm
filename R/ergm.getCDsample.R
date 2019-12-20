#  File R/ergm.getCDsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

ergm_CD_sample <- function(state, control, theta=NULL, 
                           verbose=FALSE,..., eta=ergm.eta(theta, (if(is.ergm_state(state))as.ergm_model(state)else as.ergm_model(state[[1]]))$etamap)){
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.ergm_state(state)) state <- list(state)
  states <- rep(state, length.out=nthreads(control))

  control.parallel <- control
  control.parallel$CD.samplesize <- NVL3(control$CD.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  doruns <- function(samplesize=NULL){
    if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control), ergm_CD_slave,
                                  state=states, MoreArgs=list(eta=eta,control=control.parallel,verbose=verbose,...,samplesize=samplesize))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_CD_slave(state=states[[1]], samplesize=samplesize,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  outl <- doruns()
  for(i in seq_along(outl)){
    outl[[i]]$s <- mcmc(outl[[i]]$s)
  }
  
  if(control.parallel$MCMC.runtime.traceplot){
      esteq <- lapply.mcmc.list(lapply(outl, function(out)
                      NVL3(theta, ergm.estfun(out$s, ., as.ergm_model(states[[1]])), out$s[,!as.ergm_model(states[[1]])$offsetmap,drop=FALSE])
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

ergm_CD_slave <- function(state, eta,control,verbose,..., samplesize=NULL){  
  if(is.null(samplesize)) samplesize <- control$CD.samplesize
  
  z <- if(!is.valued(state))
         .Call("CD_wrapper",
               state,
               # MCMC settings
               as.double(deInf(eta)),
               as.integer(samplesize),
               as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
               as.integer(verbose),
               PACKAGE="ergm")
       else
         .Call("WtCD_wrapper",
               state,
               # MCMC settings
               as.double(deInf(eta)),
               as.integer(samplesize),
               as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
               as.integer(verbose),
               PACKAGE="ergm")

  z$s <- matrix(z$s, ncol=nparam(state,canonical=TRUE), byrow = TRUE)
  colnames(z$s) <- param_names(state, canonical=TRUE)

  z
}
