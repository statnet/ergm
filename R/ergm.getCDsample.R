#  File R/ergm.getCDsample.R in package ergm, part of the Statnet suite
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

ergm.getCDsample <- function(nw, model, MHproposal, eta0, control, 
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

  doruns <- function(prev.runs=rep(list(NULL),nthreads), burnin=NULL, samplesize=NULL, interval=NULL){
    if(!is.null(cl)) clusterMap(cl,ergm.cdslave,
                                  Clist=Clists, prev.run=prev.runs, MoreArgs=list(MHproposal=MHproposal,eta0=eta0,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval))
    else list(ergm.cdslave(Clist=Clists[[1]], prev.run=prev.runs[[1]],burnin=burnin,samplesize=samplesize,interval=interval,MHproposal=MHproposal,eta0=eta0,control=control.parallel,verbose=verbose,...))
  }
  
  outl <- doruns()
  for(i in seq_along(outl)){
    outl[[i]]$s <- mcmc(outl[[i]]$s, control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
  }
  
  if(control.parallel$MCMC.runtime.traceplot){
    esteq <- lapply(outl, function(out)
      if(all(c("theta","etamap") %in% names(list(...)))) .ergm.esteq(list(...)$theta, list(etamap=list(...)$etamap), out$s)
      else out$s[,Clists[[1]]$diagnosable,drop=FALSE]
                    )
    for (i in seq_along(esteq)) colnames(esteq[[i]]) <- names(list(...)$theta)
    plot(as.mcmc.list(lapply(lapply(esteq, mcmc), window, thin=max(1,floor(nrow(esteq)/1000))))
        ,ask=FALSE,smooth=TRUE,density=FALSE)
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
    
    final.interval <- c(final.interval, z$final.interval)
  }
  
  ergm.stopCluster(cl)

  statsmatrix <- do.call(rbind,statsmatrices)
  colnames(statsmatrix) <- model$coef.names

  if(verbose){message("Sample size = ",nrow(statsmatrix)," by ",
                  control.parallel$MCMC.samplesize,".")}
  
  statsmatrix[is.na(statsmatrix)] <- 0
  list(statsmatrix=statsmatrix, statsmatrices=statsmatrices, status=0)

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
#
# --RETURNED--
#   the MCMC sample as a list of the following:
#     s         : the statsmatrix
#     newnwtails: the vector of tails for the new network- is this the final
#                 network sampled? - is this the original nw if 'maxedges' is 0
#     newnwheads: the vector of heads for the new network - same q's
#
###############################################################################

ergm.cdslave <- function(Clist,MHproposal,eta0,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL) {

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

  samplesize <- control$MCMC.samplesize
  
  if(is.null(Clist$weights)){
    z <- .C("CD_wrapper",
            as.integer(numnetworks), as.integer(nedges),
            as.integer(tails), as.integer(heads),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$pkgname),
            as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
            as.integer(samplesize), as.integer(c(control$CD.nsteps,control$CD.multiplicity)),
            s = as.double(rep(stats, samplesize)),
            as.integer(verbose), as.integer(MHproposal$arguments$constraints$bd$attribs),
            as.integer(MHproposal$arguments$constraints$bd$maxout), as.integer(MHproposal$arguments$constraints$bd$maxin),
            as.integer(MHproposal$arguments$constraints$bd$minout), as.integer(MHproposal$arguments$constraints$bd$minin),
            as.integer(MHproposal$arguments$constraints$bd$condAllDegExact), as.integer(length(MHproposal$arguments$constraints$bd$attribs)),
            status = integer(1),
            PACKAGE="ergm")
    
    # save the results
    z<-list(s=matrix(z$s, ncol=Clist$nstats, byrow = TRUE), status=z$status)
  }else{
    z <- .C("WtCD_wrapper",
            as.integer(length(nedges)), as.integer(nedges),
            as.integer(tails), as.integer(heads), as.double(weights),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms),
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(MHproposal$name), as.character(MHproposal$pkgname),
            as.double(c(Clist$inputs,MHproposal$inputs)), as.double(.deinf(eta0)),
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
