#  File R/ergm.getMCMCsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

#' Internal Function to Sample Networks and Network Statistics
#' 
#' This is an internal function, not normally called directly by the
#' user. The \code{ergm_MCMC_sample} function samples networks and
#' network statistics using an MCMC algorithm via \code{MCMC_wrapper}
#' and is caple of running in multiple threads using
#' `ergm_MCMC_slave`.
#' 
#' 
#' Note that unless `stats0` is passed, the returned stats will be relative to the original network, i.e.,
#' the calling function must shift the statistics if required.
#' 
#' @param nw a [`network`] (or [`pending_update_network`]) object representing the sampler state.
#' @param model an [`ergm_model`] to be sampled from, as returned by
#'   [ergm_model()].
#' @param proposal a list of the parameters needed for
#'   Metropolis-Hastings proposals and the result of calling
#'   [ergm_proposal()].
#' @param control list of MCMC tuning parameters; see
#'   [control.ergm()].
#' @param theta the (possibly curved) parameters of the model.
#' @param eta the natural parameters of the model; by default constructed from `theta`.
#' @param stats0 either a numeric vector of the same length as `eta`
#'   containing the statistics corresponding to the initial network;
#'   or a list thereof of the same length as that of `nw` if it is a
#'   list of networks (for parallel sampling). Defaults to a vector of
#'   0s. Note that if passed to `ergm_MCMC_slave`, it is overridden if
#'   `prev.run` is passed as well.
#' @param verbose verbosity level.
#' @template response
#' @param update.nws whether to actually update the network state or
#'   to return an object "promising" to update the network.
#' @param ... additional arugments.
#'
#' @return
#' \code{ergm_MCMC_sample} returns a list
#'   containing:
#' \item{stats}{an [`mcmc.list`] with sampled statistics.}
#' \item{networks}{a list of final sampled networks, one for each thread.}
#' \item{status}{status code, propagated from `ergm.mcmcslave`.}
#' \item{final.interval}{adaptively determined MCMC interval.}
#'
#' If `update.nws==FALSE`, rather than returning the updated networks,
#' the function will return a [`pending_update_network`].
#'
#' @note `ergm_MCMC_sample` and `ergm_MCMC_slave` replace
#'   `ergm.getMCMCsample` and `ergm.mcmcslave` respectively. They
#'   differ slightly in their argument names and in their return
#'   formats. For example, `ergm_MCMC_sample` expects `proposal`
#'   rather than `MHproposal` and `theta` or `eta` rather than `eta0`;
#'   and it does not return `statsmatrix` or `newnetwork`
#'   elements. Rather, if parallel processing is not in effect,
#'   `stats` is an [`mcmc.list`] with one chain and `networks` is a
#'   list with one element.
#' @export
ergm_MCMC_sample <- function(nw, model, proposal, control, theta=NULL, 
                             response=NULL, update.nws = TRUE, verbose=FALSE,..., eta=ergm.eta(theta, model$etamap), stats0=NULL) {
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.network(nw) || is.pending_update_network(nw)) nw <- list(nw)
  nws <- rep(nw, length.out=nthreads(control))
  NVL(stats0) <- numeric(length(eta))
  if(is.numeric(stats0)) stats0 <- list(stats0)
  stats0 <- rep(stats0, length.out=length(nws))

  Clists <- lapply(nws, ergm::ergm.Cprepare, model, response=response)

  control.parallel <- control
  control.parallel$MCMC.samplesize <- NVL3(control$MCMC.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  #' @importFrom parallel clusterMap
  doruns <- function(prev.runs=rep(list(NULL),nthreads(control)), burnin=NULL, samplesize=NULL, interval=NULL){
    if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control),ergm_MCMC_slave,
                                  Clist=Clists, prev.run=prev.runs, stats0=stats0, MoreArgs=list(proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_MCMC_slave(Clist=Clists[[1]], prev.run=prev.runs[[1]], stats0=stats0[[1]], burnin=burnin,samplesize=samplesize,interval=interval,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  if(!is.null(control.parallel$MCMC.effectiveSize)){
    if(verbose) message("Beginning adaptive MCMC...")

    howmuchmore <- function(target.ess, current.ss, current.ess, current.burnin){
      (target.ess-current.ess)*(current.ss-current.burnin)/current.ess
    }
    
    interval <- control.parallel$MCMC.interval
    best.burnin <- list(burnin=0,pval=control.parallel$MCMC.effectiveSize.burnin.pval)
    eS <- 0
    outl <- rep(list(NULL),nthreads(control))
    for(mcrun in seq_len(control.parallel$MCMC.effectiveSize.maxruns)){
      if(mcrun==1){
        samplesize <- control.parallel$MCMC.samplesize
        if(verbose>1)
          message("First run: running each chain forward by ",samplesize, " steps with interval ", interval, ".")
      }else{
        if(burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
          samplesize <- control.parallel$MCMC.samplesize
          if(verbose>1)
            message("Insufficient ESS or untrustworthy burn-in estimate to determine the number of steps remaining: running forward by ",samplesize, " steps with interval ", interval, ".")
        }else{
          pred.ss <- howmuchmore(control.parallel$MCMC.effectiveSize, NVL(nrow(outl[[1]]$s),0), eS, best.burnin$burnin)
          damp.ss <- pred.ss*(eS/(control.parallel$MCMC.effectiveSize.damp+eS))+control.parallel$MCMC.samplesize*(1-eS/(control.parallel$MCMC.effectiveSize.damp+eS))
          samplesize <- round(damp.ss)
          if(verbose>1) message("Predicted additional sample size: ",pred.ss, " dampened to ",damp.ss, ", so running ", samplesize, " steps forward.")
        }
      }
        
      outl<-doruns(prev.runs=outl,
                  burnin = interval, # I.e., skip that much before the first draw.
                  samplesize = samplesize,
                  interval = interval
                  )
      
      # Stop if something went wrong.
      if(any(sapply(outl,"[[","status")!=0)) break
      
      while(nrow(outl[[1]]$s)-best.burnin$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)) outl[[i]]$s <- outl[[i]]$s[seq_len(floor(nrow(outl[[i]]$s)/2))*2+nrow(outl[[i]]$s)%%2,,drop=FALSE]
        interval <- interval*2
        if(verbose) message("Increasing thinning to ",interval,".")
      }
      
      esteq <- lapply.mcmc.list(lapply(outl, function(out)
                      NVL3(theta, ergm.estfun(out$s, ., model), out$s[,Clists[[1]]$diagnosable,drop=FALSE])
                      ), mcmc, start=1, thin=interval)
      
      if(control.parallel$MCMC.runtime.traceplot){
        plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
             ,ask=FALSE,smooth=TRUE,density=FALSE)
      }

      best.burnin <- .find_OK_burnin(esteq, npts=control$MCMC.effectiveSize.points, base=control$MCMC.effectiveSize.base, min.pval=control$MCMC.effectiveSize.burnin.pval, order.max=control$MCMC.effectiveSize.order.max)
      burnin.pval <- best.burnin$pval
      if(burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
        if(verbose>1) message("No adequate burn-in found. Best convergence p-value = ", burnin.pval)
        next
      }
      postburnin.mcmc <- window(esteq, start=start(esteq)+best.burnin$burnin*thin(esteq))
      
      eS <- niter(postburnin.mcmc)*nchain(postburnin.mcmc)/attr(spectrum0.mvar(postburnin.mcmc, order.max=control$MCMC.effectiveSize.order.max),"infl")
      
      if(verbose) message("ESS of ",eS," attained with burn-in of ", round(best.burnin$burnin/niter(esteq)*100,2),"%; convergence p-value = ", burnin.pval, ".")

      if(eS>=control.parallel$MCMC.effectiveSize){
        if(burnin.pval > control$MCMC.effectiveSize.burnin.pval){
          if(verbose) message("Target ESS achieved and is trustworthy. Returning.")
          break
        }else{
          if(verbose>1) message("ESS and burn-in estimates are not trustworthy.")
        }
      }
    }

    if(eS<control.parallel$MCMC.effectiveSize)
      warning("Unable to reach target effective size in iterations alotted.")

    for(i in seq_along(outl)){
      if(best.burnin$burnin) outl[[i]]$s <- outl[[i]]$s[-seq_len(best.burnin$burnin),,drop=FALSE]
      outl[[i]]$s <- coda::mcmc(outl[[i]]$s, (best.burnin$burnin+1)*interval, thin=interval)
      outl[[i]]$final.interval <- interval
    }
  }else{
    outl <- doruns()
    for(i in seq_along(outl)){
      outl[[i]]$s <- coda::mcmc(outl[[i]]$s, control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
    }
    
    if(control.parallel$MCMC.runtime.traceplot){
      esteq <- as.mcmc.list(lapply(lapply(outl, function(out)
        NVL3(theta, ergm.estfun(out$s, ., model), out$s[,Clists[[1]]$diagnosable,drop=FALSE])
      ), mcmc, start=control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval))
      plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
           ,ask=FALSE,smooth=TRUE,density=FALSE)
    }
  }

  #
  #   Process the results
  #
  statsmatrices <- list()
  newnetworks <- list()
  final.interval <- c()
  for(i in (1:nthreads(control))){
    z <- outl[[i]]
    
    if(z$status == 1){ # MCMC_TOO_MANY_EDGES, exceeding even control.parallel$MCMC.maxedges
      return(list(status=z$status))
    }
    
    if(z$status == 2){ # MCMC_MH_FAILED
      # MH proposal failed somewhere. Throw an error.
      stop("Sampling failed due to a Metropolis-Hastings proposal failing.")
      }
    
    statsmatrices[[i]] <- z$s

    newnetworks[[i]] <- pending_update_network(nws[[i]], z, response=response)
    if(update.nws){
      newnetworks[[i]] <- as.network(newnetworks[[i]])
    }
    final.interval <- c(final.interval, z$final.interval)
  }
  
  stats <- as.mcmc.list(statsmatrices)
  if(verbose){message("Sample size = ",niter(stats)*nchain(stats)," by ",
                  niter(stats),".")}
  
  list(stats = stats, networks=newnetworks, status=0, final.interval=final.interval)
}

#' @rdname ergm_MCMC_sample
#' @description The \code{ergm_MCMC_slave} function calls the actual C
#'   routine and does minimal preprocessing.
#'
#' @param prev.run a summary of the state of the sampler allowing a
#'   run to be resumed quickly by `ergm_MCMC_slave`.
#' @param burnin,samplesize,interval MCMC paramters that can
#'   be used to temporarily override those in the `control` list.
#' @param Clist the list of parameters returned by
#'   \code{\link{ergm.Cprepare}}
#' @return \code{ergm_MCMC_slave} returns the MCMC sample as a list of
#'   the following: \item{s}{the matrix of statistics.}
#'   \item{newnwtails}{the vector of tails for the new network.}
#'   \item{newnwheads}{the vector of heads for the new network.}
#'   \item{newnwweights}{the vector of weights for the new network (if
#'   applicable)} \item{status}{success or failure code: `0` is
#'   success, `1` for too many edges, and `2` for a
#'   Metropolis-Hastings proposal failing.}
#' @useDynLib ergm
#' @export
ergm_MCMC_slave <- function(Clist,proposal,eta,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, stats0=numeric(Clist$nstats)) {
  if(is.null(prev.run)){ # Start from Clist
    nedges <- Clist$nedges
    tails <- Clist$tails
    heads <- Clist$heads
    weights <- Clist$weights
    stats <- stats0
  }else{ # Pick up where we left off
    nedges <- prev.run$newnwtails[1]
    tails <- prev.run$newnwtails[2:(nedges+1)]
    heads <- prev.run$newnwheads[2:(nedges+1)]
    weights <- prev.run$newnwweights[2:(nedges+1)]
    stats <- prev.run$s[nrow(prev.run$s),]
  }

  if(is.null(burnin)) burnin <- control$MCMC.burnin
  if(is.null(samplesize)) samplesize <- control$MCMC.samplesize
  if(is.null(interval)) interval <- control$MCMC.interval

  z <-
    if(is.null(Clist$weights)){
      .Call("MCMC_wrapper",
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
            as.integer(burnin), 
            as.integer(interval),
            as.integer(.deinf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")
    }else{
      .Call("WtMCMC_wrapper",
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
            as.integer(burnin), 
            as.integer(interval),
            as.integer(.deinf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")
    }

  # save the results (note that if prev.run is NULL, c(NULL$s,z$s)==z$s.
  z$s <- matrix(z$s+stats, ncol=Clist$nstats, byrow = TRUE)
  z$s <- rbind(prev.run$s,z$s)
  colnames(z$s) <- names(Clist$diagnosable)

  z
}


.find_OK_burnin <- function(x, npts, base, min.pval=0.2, ...){
  geweke <- function(b){
    if(b>0) x <- window(x, start=start(x) + b * thin(x))
    suppressWarnings(geweke.diag.mv(x, ...)$p.value)
  }

  # TODO: Implement bisection algorithm here.
  pts <- sort(round(base^seq_len(npts)*niter(x)))
  pvals <- sapply(pts, geweke)

  best <- suppressWarnings(min(which(pvals>min.pval)))
  if(!is.finite(best))
    best <- which.max(pvals)
  
  list(burnin=pts[best], pval=pvals[best])
}
