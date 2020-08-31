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
#' @param state an [`ergm_state`] representing the sampler state, containing information about the network, the model, the proposal, and (optionally) initial statistics, or a list thereof.
#' @param control list of MCMC tuning parameters; see
#'   [control.ergm()].
#' @param theta the (possibly curved) parameters of the model.
#' @param eta the natural parameters of the model; by default constructed from `theta`.
#' @param verbose verbosity level.
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
#' the function will return a [`ergm_state`].
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
ergm_MCMC_sample <- function(state, control, theta=NULL, 
                             verbose=FALSE,..., eta=ergm.eta(theta, (if(is.ergm_state(state))as.ergm_model(state)else as.ergm_model(state[[1]]))$etamap)){
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.ergm_state(state)) state <- list(state)
  state <- rep(state, length.out=nthreads(control))
  ## state <- if(is.network(state[[1]])){
  ##   NVL(stats0) <- numeric(length(eta))
  ##   if(is.numeric(stats0)) stats0 <- list(stats0)
  ##   stats0 <- rep(stats0, length.out=length(state))

  ##   state <- mapply(ergm_state, state, stats=stats0, MoreArgs=list(response=response, model=model, proposal=proposal), SIMPLIFY=FALSE)
  ## }else state

  control.parallel <- control
  control.parallel$MCMC.samplesize <- NVL3(control$MCMC.samplesize, ceiling(. / nthreads(control)))

  flush.console()

  state0 <- state
  state <- lapply(state, ergm_state_send) # Don't carry around nw0.
  
  #' @importFrom parallel clusterMap
  doruns <- function(burnin=NULL, samplesize=NULL, interval=NULL){
    out <- if(!is.null(ergm.getCluster(control))) persistEvalQ({clusterMap(ergm.getCluster(control),ergm_MCMC_slave,
                                                                    state=state, MoreArgs=list(eta=eta,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval))}, retries=getOption("ergm.cluster.retries"), beforeRetry={ergm.restartCluster(control,verbose)})
    else list(ergm_MCMC_slave(state[[1]], burnin=burnin,samplesize=samplesize,interval=interval,eta=eta,control=control.parallel,verbose=verbose,...))
    # Note: the return value's state will be a ergm_state_receive.
    for(i in seq_along(out)) out[[i]]$state <- update(state[[i]], out[[i]]$state)
    out
  }

  sms <- vector("list", nthreads(control))
  
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
          pred.ss <- howmuchmore(control.parallel$MCMC.effectiveSize, NVL(nrow(sms[[1]]),0), eS, best.burnin$burnin)
          damp.ss <- pred.ss*(eS/(control.parallel$MCMC.effectiveSize.damp+eS))+control.parallel$MCMC.samplesize*(1-eS/(control.parallel$MCMC.effectiveSize.damp+eS))
          samplesize <- round(damp.ss)
          if(verbose>1) message("Predicted additional sample size: ",pred.ss, " dampened to ",damp.ss, ", so running ", samplesize, " steps forward.")
        }
      }
        
      outl<-doruns(burnin = interval, # I.e., skip that much before the first draw.
                   samplesize = samplesize,
                   interval = interval)

      # Stop if something went wrong.
      if(any(map_int(outl,"status")!=0)) break
      sms <- mapply(rbind, sms, map(outl, "s"), SIMPLIFY=FALSE)
      state <- map(outl, "state")
      
      while(nrow(sms[[1]])-best.burnin$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)) sms[[i]] <- sms[[i]][seq_len(floor(nrow(sms[[i]])/2))*2+nrow(sms[[i]])%%2,,drop=FALSE]
        interval <- interval*2
        if(verbose) message("Increasing thinning to ",interval,".")
      }
      
      esteq <- lapply(sms, function(sm) NVL3(theta, ergm.estfun(sm, ., as.ergm_model(state[[1]])), sm[,!as.ergm_model(state[[1]])$etamap$offsetmap,drop=FALSE])) %>%
        lapply.mcmc.list(mcmc, start=1, thin=interval) %>% lapply.mcmc.list(`-`)

      if(control.parallel$MCMC.runtime.traceplot){
        plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
            ,ask=FALSE,smooth=TRUE,density=FALSE)
      }

      best.burnin <- .find_OK_burnin(esteq, order.max=control$MCMC.effectiveSize.order.max)
      burnin.pval <- best.burnin$pval
      if(burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
        if(verbose>1) message("Selected burn-in p-value = ", burnin.pval, " is below the threshold of ",control$MCMC.effectiveSize.burnin.pval,".")
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
      if(best.burnin$burnin) sms[[i]] <- sms[[i]][-seq_len(best.burnin$burnin),,drop=FALSE]
      sms[[i]] <- coda::mcmc(sms[[i]], (best.burnin$burnin+1)*interval, thin=interval)
      outl[[i]]$final.interval <- interval
    }
  }else{
    outl <- doruns()
    sms <- map(outl, "s")
    for(i in seq_along(outl)){
      sms[[i]] <- coda::mcmc(sms[[i]], control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
    }
    
    if(control.parallel$MCMC.runtime.traceplot){
      lapply(sms, function(sm) NVL3(theta, ergm.estfun(sm, ., as.ergm_model(state[[1]])), sm[,!as.ergm_model(state[[1]])$etamap$offsetmap,drop=FALSE])) %>% lapply(mcmc, start=control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval) %>% as.mcmc.list() %>% window(., thin=thin(.)*max(1,floor(niter(.)/1000))) %>% plot(ask=FALSE,smooth=TRUE,density=FALSE)
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
    
    statsmatrices[[i]] <- sms[[i]]

    newnetworks[[i]] <- update(state0[[i]], state=z$state)
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
#' @param burnin,samplesize,interval MCMC paramters that can
#'   be used to temporarily override those in the `control` list.
#' @return \code{ergm_MCMC_slave} returns the MCMC sample as a list of
#'   the following: \item{s}{the matrix of statistics.}
#'   \item{state}{an [`ergm_state`] object for the new network.}
#'   \item{status}{success or failure code: `0` is
#'   success, `1` for too many edges, and `2` for a
#'   Metropolis-Hastings proposal failing.}
#' @useDynLib ergm
#' @export
ergm_MCMC_slave <- function(state, eta,control,verbose,..., burnin=NULL, samplesize=NULL, interval=NULL){
  on.exit(ergm_Cstate_clear())

  NVL(burnin) <- control$MCMC.burnin
  NVL(samplesize) <- control$MCMC.samplesize
  NVL(interval) <- control$MCMC.interval
  z <-
    if(!is.valued(state))
      .Call("MCMC_wrapper",
            state,
            # MCMC settings
            as.double(deInf(eta)),
            as.integer(samplesize),
            as.integer(burnin), 
            as.integer(interval),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")
    else
      .Call("WtMCMC_wrapper",
            state,
            # MCMC settings
            as.double(deInf(eta)),
            as.integer(samplesize),
            as.integer(burnin), 
            as.integer(interval),
            as.integer(deInf(NVL(control$MCMC.maxedges,Inf),"maxint")),
            as.integer(verbose),
            PACKAGE="ergm")

  if(z$status) return(z) # If there is an error.
  z$s <- matrix(z$s, ncol=nparam(state,canonical=TRUE), byrow = TRUE)
  colnames(z$s) <- param_names(state, canonical=TRUE)
  z$state <- ergm_state_receive(z$state)

  z
}


.find_OK_burnin <- function(x, ...){
  xs <- x %>% map(scale) %>% map(~.[,attr(.,"scaled:scale")>0,drop=FALSE])
  ssr <- function(b,s){
    b <- round(b)
    n <- nrow(s)
    sum(resid(lm(s~c(seq_len(b), numeric(n-b))+ I(seq_len(n)>b)))^2)
  }

  geweke <- function(b){
    if(b>0) x <- window(x, start=start(x) + b * thin(x))
    p.val <- suppressWarnings(geweke.diag.mv(x, ...)$p.value)
    if(is.na(p.val)) 0 else p.val
  }

  best <- max(sapply(xs, function(x) optimize(ssr, c(0, nrow(x)/2), s=x, tol=1)$minimum))

  list(burnin=round(best), pval=geweke(round(best)))
}
