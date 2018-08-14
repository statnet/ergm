#  File R/ergm.getMCMCsample.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################

#' @include ergm-deprecated.R
#' @describeIn ergm-deprecated Use [ergm_MCMC_sample()] instead.
#' @export ergm.getMCMCsample
ergm.getMCMCsample <- function(nw, model, proposal, eta0, control, 
                               verbose=FALSE, response=NULL, update.nws = TRUE,...) {
  .dep_once("ergm_MCMC_sample")
  out <- ergm_MCMC_sample(nw, model, proposal, eta=eta0, control=control, verbose=verbose, response=response, update.nws=update.nws, theta=list(...)$theta, ...)
  
  out$newnetworks<-out$networks
  out$newnetwork<-out$newnetworks[[1]]
  out$statsmatrices <- out$stats
  out$stats <- NULL
  out$statsmatrix <- do.call(rbind,out$statsmatrices)
  colnames(out$statsmatrix) <- model$coef.names
  out$statsmatrix[is.na(out$statsmatrix)] <- 0
  out
}

#' Internal Function to Sample Networks and Network Statistics
#' 
#' This is an internal function, not normally called directly by the
#' user. The \code{ergm_MCMC_sample} function samples networks and
#' network statistics using an MCMC algorithm via \code{MCMC_wrapper}
#' and is caple of running in multiple threads using
#' `ergm_MCMC_slave`.
#' 
#' 
#' Note that the returned stats will be relative to the original network, i.e.,
#' the calling function must shift the statistics if required. The calling
#' function must also attach column names to the statistics matrix if required.
#' 
#' @param nw a [`network`] object representing the sampler state.
#' @param model an [`ergm_model`] to be sampled from, as returned by
#'   [ergm_model()].
#' @param proposal a list of the parameters needed for
#'   Metropolis-Hastings proposals and the result of calling
#'   [ergm_proposal()].
#' @param control list of MCMC tuning parameters; see
#'   [control.ergm()].
#' @param theta the (possibly curved) parameters of the model.
#' @param eta the natural parameters of the model; by default constructed from `theta`.
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
#' the function will remove all edges from the input networks, attach
#' a network attribute `.update` with the new edge information, and
#' change class name to prevent the resulting object from being
#' accessed or modified by functions that do not understand it.
#'
#' @note Unlike its predecessor `ergm.getMCMCsample`,
#'   `ergm_MCMC_sample` does not return `statsmatrix` or `newnetwork`
#'   elements. Rather, if parallel processing is not in effect,
#'   `stats` is an [`mcmc.list`] with one chain and
#'   `networks` is a list with one element.
#' @export ergm_MCMC_sample
ergm_MCMC_sample <- function(nw, model, proposal, control, theta=NULL, 
                             response=NULL, update.nws = TRUE, verbose=FALSE,..., eta=ergm.eta(theta, model$etamap)) {
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
  control.parallel$MCMC.samplesize <- NVL3(control$MCMC.samplesize, ceiling(. / nthreads))

  flush.console()

  #' @importFrom parallel clusterMap
  doruns <- function(prev.runs=rep(list(NULL),nthreads), burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL){
    if(!is.null(cl)) clusterMap(cl,ergm_MCMC_slave,
                                  Clist=Clists, prev.run=prev.runs, MoreArgs=list(proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval,maxedges=maxedges))
    else list(ergm_MCMC_slave(Clist=Clists[[1]], prev.run=prev.runs[[1]],burnin=burnin,samplesize=samplesize,interval=interval,maxedges=maxedges,proposal=proposal,eta=eta,control=control.parallel,verbose=verbose,...))
  }
  
  if(!is.null(control.parallel$MCMC.effectiveSize)){
    if(verbose) message("Beginning adaptive MCMC...")

    howmuchmore <- function(target.ess, current.ss, current.ess, current.burnin){
      (target.ess-current.ess)*(current.ss-current.burnin)/current.ess
    }
    
    interval <- control.parallel$MCMC.interval
    meS <- list(burnin=0,eS=control.parallel$MCMC.effectiveSize)
    outl <- rep(list(NULL),nthreads)
    for(mcrun in seq_len(control.parallel$MCMC.effectiveSize.maxruns)){
      if(mcrun==1){
        samplesize <- control.parallel$MCMC.samplesize
        if(verbose>1)
          message("First run: running each chain forward by ",samplesize, " steps with interval ", interval, ".")
      }else{
        if(meS$eS<1 || meS$pts.rank==1){
          samplesize <- control.parallel$MCMC.samplesize
          if(verbose>1)
            message("Insufficient ESS or untrustworthy burn-in estimate to determine the number of steps remaining: running forward by ",samplesize, " steps with interval ", interval, ".")
        }else{
          pred.ss <- howmuchmore(control.parallel$MCMC.effectiveSize, NVL(nrow(outl[[1]]$s),0), meS$eS, meS$burnin)
          damp.ss <- pred.ss*(meS$eS/(control.parallel$MCMC.effectiveSize.damp+meS$eS))+control.parallel$MCMC.samplesize*(1-meS$eS/(control.parallel$MCMC.effectiveSize.damp+meS$eS))
          samplesize <- round(damp.ss)
          if(verbose>1) message("Predicted additional sample size: ",pred.ss, " dampened to ",damp.ss, ", so running ", samplesize, " steps forward.")
        }
      }
        
      outl<-doruns(prev.runs=outl,
                  burnin = interval, # I.e., skip that much before the first draw.
                  samplesize = samplesize,
                  interval = interval,
                  maxedges = if(!is.null(outl[[1]])) max(sapply(outl,"[[","maxedges"))
                  )
      
      # Stop if something went wrong.
      if(any(sapply(outl,"[[","status")!=0)) break
      
      while(nrow(outl[[1]]$s)-meS$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)) outl[[i]]$s <- outl[[i]]$s[seq_len(floor(nrow(outl[[i]]$s)/2))*2+nrow(outl[[i]]$s)%%2,,drop=FALSE]
        interval <- interval*2
        if(verbose) message("Increasing thinning to ",interval,".")
      }
      
      esteq <- as.mcmc.list(lapply(lapply(outl, function(out)
                      NVL3(theta, ergm.estfun(out$s, ., model), out$s[,Clists[[1]]$diagnosable,drop=FALSE])
                      ), mcmc, start=1, thin=interval))
      
      meS <- .max.effectiveSize(esteq, npts=control$MCMC.effectiveSize.points, base=control$MCMC.effectiveSize.base, ar.order=control$MCMC.effectiveSize.order)
      if(verbose>1) message("Maximum harmonic mean ESS of ",meS$eS," attained with burn-in of ", round(meS$b/nrow(outl[[1]]$s)*100,2),"%.")

      if(control.parallel$MCMC.runtime.traceplot){
        plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
             ,ask=FALSE,smooth=TRUE,density=FALSE)
      }

      if(meS$pts.rank!=1 && meS$eS>=control.parallel$MCMC.effectiveSize){
        if(verbose) message("Target ESS achieved and is trustworthy. Returning.")
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

    if(update.nws){
      newnetworks[[i]] <- newnw.extract(nws[[i]],z, response=response,
                                        output=control.parallel$network.output)
    }else{
      if(length(newnetworks)<i) newnetworks[[i]] <- nws[[i]]
      class(newnetworks[[i]]) <- "network"
      if(network.edgecount(newnetworks[[i]])!=0) newnetworks[[i]] <- empty_network(newnetworks[[i]])
      newnetworks[[i]]%n%".update" <- z[c("newnwtails","newnwheads","newnwweights")]
      # Make sure that nobody treats this as an actual network object
      # by mistake.
      class(newnetworks[[i]]) <- "pending_update_network"
    }
    final.interval <- c(final.interval, z$final.interval)
  }
  
  ergm.stopCluster(cl)

  stats <- as.mcmc.list(statsmatrices)
  if(verbose){message("Sample size = ",niter(stats)*nchain(stats)," by ",
                  niter(stats),".")}
  
  list(stats = stats, networks=newnetworks, status=0, final.interval=final.interval)
}

#' @include ergm-deprecated.R
#' @describeIn ergm-deprecated Use [ergm_MCMC_slave()] instead.
#' @export ergm.mcmcslave
ergm.mcmcslave <- function(Clist,proposal,eta0,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL){
  .dep_once("ergm_MCMC_slave")
  ergm_MCMC_slave(Clist,proposal,eta0,control,verbose,...,prev.run=prev.run, burnin=burnin, samplesize=samplesize, interval=interval, maxedges=maxedges)
}


#' @rdname ergm_MCMC_sample
#' @description The \code{ergm_MCMC_slave} function calls the actual C
#'   routine and does minimal preprocessing.
#'
#' @param prev.run a summary of the state of the sampler allowing a
#'   run to be resumed quickly by `ergm_MCMC_slave`.
#' @param burnin,samplesize,interval,maxedges MCMC paramters that can
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
#'   Metropolis-Hastings proposal failing.}  \item{maxedges}{maximum
#'   allowed edges at the time of return.}
#' @useDynLib ergm
#' @export ergm_MCMC_slave
ergm_MCMC_slave <- function(Clist,proposal,eta,control,verbose,...,prev.run=NULL, burnin=NULL, samplesize=NULL, interval=NULL, maxedges=NULL) {
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
              as.character(proposal$name), as.character(proposal$pkgname),
              as.double(c(Clist$inputs,proposal$inputs)), as.double(.deinf(eta)),
              as.integer(samplesize),
              s = as.double(rep(stats, samplesize)),
              as.integer(burnin), 
              as.integer(interval),
              newnwtails = integer(maxedges),
              newnwheads = integer(maxedges),
              as.integer(verbose), as.integer(proposal$arguments$constraints$bd$attribs),
              as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
              as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
              as.integer(proposal$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal$arguments$constraints$bd$attribs)),
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
              as.character(proposal$name), as.character(proposal$pkgname),
              as.double(c(Clist$inputs,proposal$inputs)), as.double(.deinf(eta)),
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
    if(ar.order) .fast.effectiveSize(as.matrix(coda::as.mcmc.list(x), ar.order=ar.order))
    else effectiveSize(as.matrix(coda::as.mcmc.list(x)))
  }

  # TODO: Implement bisection algorithm here.
  pts <- sort(round(base^seq_len(npts)*nrow(x[[1]])))
  ess <- rbind(sapply(pts, es)) # I.e., variables in rows and burn-in test points in columns.

  best <- max(apply(ess, 1, which.max))

  mean.fn <- function(x) x^(-1)
  mean.ifn <- function(x) x^(-1)
  hmean <- mean.ifn(mean(mean.fn(ess[,best])))
  
  list(burnin=pts[best], eS=hmean, pts.rank=length(pts)-best+1)
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
      novar <- var(x[,i])<.Machine$double.eps
      if(novar){
        v0[i] <- 0
        order[i] <- 0
      }else{
        ar.out <- ar(x[, i], aic = FALSE, order.max=ar.order)
        v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
        order[i] <- ar.out$order
      }
    }
    return(list(spec = v0, order = order))
}
