#  File R/ergm.getCDsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

ergm_CD_sample <- function(nw, model, proposal, control, theta=NULL, 
                             response=NULL, verbose=FALSE,..., eta=ergm.eta(theta, model$etamap)) {
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.network(nw) || is.pending_update_network(nw)) nw <- list(nw)
  nws <- rep(nw, length.out=nthreads(control))
  
  Clists <- lapply(nws, ergm::ergm.Cprepare, model, response=response)

  control.parallel <- control
  control.parallel$MCMC.samplesize <- NVL3(control$MCMC.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  doruns <- function(samplesize=NULL){
    if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control), ergm_CD_slave,
                                  Clist=Clists, MoreArgs=list(proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...,samplesize=samplesize))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_CD_slave(Clist=Clists[[1]], samplesize=samplesize,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  outl <- doruns()
  for(i in seq_along(outl)){
    outl[[i]]$s <- mcmc(outl[[i]]$s)
  }
  
  if(control.parallel$MCMC.runtime.traceplot){
    lapply(outl, function(out) NVL3(theta, ergm.estfun(out$s, ., model), out$s[,Clists[[1]]$diagnosable,drop=FALSE])) %>%
      lapply.mcmc.list(mcmc, start=1) %>% lapply.mcmc.list(`-`) %>% window(., thin=thin(.)*max(1,floor(niter(.)/1000))) %>%
      plot(ask=FALSE,smooth=TRUE,density=FALSE)
  }

  #
  #   Process the results
  #
  statsmatrices <- list()
  newnetworks <- list()
  for(i in (1:nthreads(control))){
    z <- outl[[i]]
    
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control.parallel$MCMC.max.maxedges
      return(list(status=z$status))
    }
    
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

ergm_CD_slave <- function(Clist,proposal,eta,control,verbose,..., samplesize=NULL) {
    nedges <- c(Clist$nedges,0,0)
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
    stats <- rep(0, Clist$nstats)
  
  if(is.null(samplesize)) samplesize <- control$MCMC.samplesize

  samplesize <- control$MCMC.samplesize
  
  if(is.null(Clist$weights)){
    z <- .C("CD_wrapper",
            as.integer(nedges),
            as.integer(tails), as.integer(heads),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(proposal$name), as.character(proposal$pkgname),
            as.double(c(Clist$inputs,proposal$inputs)), as.double(deInf(eta)),
            as.integer(samplesize), as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
            s = as.double(rep(stats, samplesize)),
            as.integer(verbose), as.integer(proposal$arguments$constraints$bd$attribs),
            as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
            as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
            as.integer(proposal$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal$arguments$constraints$bd$attribs)),
            status = integer(1),
            PACKAGE="ergm")
    
    # save the results
    z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE), status=z$status)
  }else{
    z <- .C("WtCD_wrapper",
            as.integer(nedges),
            as.integer(tails), as.integer(heads), as.double(weights),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(proposal$name), as.character(proposal$pkgname),
            as.double(c(Clist$inputs,proposal$inputs)), as.double(deInf(eta)),
            as.integer(samplesize), as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
            s = as.double(rep(stats, samplesize)),
            as.integer(verbose), 
            status = integer(1),
            PACKAGE="ergm")
    # save the results
    z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE), status=z$status)
  }
  
  z
    
}
