#  File R/ergm.getMCMCsample.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' Internal Function to Sample Networks and Network Statistics
#' 
#' This is an internal function, not normally called directly by the
#' user. The \code{ergm_MCMC_sample} function samples networks and
#' network statistics using an MCMC algorithm via \code{MCMC_wrapper}
#' and is capable of running in multiple threads using
#' `ergm_MCMC_slave`.
#' 
#' 
#' @param state an [`ergm_state`] representing the sampler state, containing information about the network, the model, the proposal, and (optionally) initial statistics, or a list thereof.
#'
#' @templateVar mycontrols [control.ergm()], [control.simulate.ergm()], etc.
#' @template control2
#' @template verbose
#'
#' @param theta the (possibly curved) parameters of the model.
#' @param eta the natural parameters of the model; by default constructed from `theta`.
#' @param ... additional arugments.
#'
#' @return
#' \code{ergm_MCMC_sample} returns a list
#'   containing:
#' \item{stats}{an [`mcmc.list`] with sampled statistics.}
#' \item{networks}{a list of final sampled networks, one for each thread.}
#' \item{status}{status code, propagated from `ergm_MCMC_slave()`.}
#' \item{final.interval}{adaptively determined MCMC interval.}
#'
#' \item{sampnetworks}{If `control$MCMC.save_networks` is set and is
#' `TRUE`, a list of lists of `ergm_state`s corresponding to the
#' sampled networks.}
#'
#' @note `ergm_MCMC_sample` and `ergm_MCMC_slave` replace
#'   `ergm.getMCMCsample` and `ergm.mcmcslave` respectively. They
#'   differ slightly in their argument names and in their return
#'   formats. For example, `ergm_MCMC_sample` expects `ergm_state`
#'   rather than network/model/proposal, and `theta` or `eta` rather than `eta0`;
#'   and it does not return `statsmatrix` or `newnetwork`
#'   elements. Rather, if parallel processing is not in effect,
#'   `stats` is an [`mcmc.list`] with one chain and `networks` is a
#'   list with one element.
#'
#'   Note that unless `stats` is a part of the `ergm_state`, the
#'   returned stats will be relative to the original network, i.e.,
#'   the calling function must shift the statistics if required.
#'
#'   At this time, repeated calls to `ergm_MCMC_sample` will not
#'   produce the same sequence of networks as a single long call, even
#'   with the same starting seeds. This is because the network
#'   sampling algorithms rely on the internal state of the network
#'   representation in C, which may not be reconstructed exactly the
#'   same way when "resuming". This behaviour may change in the
#'   future.
#'
#' @examples
#'
#' # This example illustrates constructing "ingredients" for calling
#' # ergm_MCMC_sample() from calls to simulate.ergm(). One can also
#' # construct an ergm_state object directly from ergm_model(),
#' # ergm_proposal(), etc., but the approach shown here is likely to
#' # be the least error-prone and the most robust to future API
#' # changes.
#' #
#' # The regular simulate() call hierarchy is
#' #
#' # simulate_formula.network(formula) ->
#' #   simulate.ergm_model(ergm_model) ->
#' #     simulate.ergm_state_full(ergm_state)
#' #
#' # They take an argument, return.args=, that will interrupt the call
#' # and have it return its arguments. We can use it to obtain
#' # low-level inputs robustly.
#'
#' data(florentine)
#' control <- control.simulate(MCMC.burnin = 2, MCMC.interval = 1)
#'
#'
#' # FYI: Obtain input for simulate.ergm_model():
#' sim.mod <- simulate(flomarriage~absdiff("wealth"), constraints=~edges,
#'                     coef = NULL, nsim=3, control=control,
#'                     return.args="ergm_model")
#' names(sim.mod)
#' str(sim.mod$object,1) # ergm_model
#'
#' # Obtain input for simulate.ergm_state_full():
#' sim.state <- simulate(flomarriage~absdiff("wealth"), constraints=~edges,
#'                       coef = NULL, nsim=3, control=control,
#'                       return.args="ergm_state")
#' names(sim.state)
#' str(sim.state$object, 1) # ergm_state
#'
#' # This control parameter would be set by nsim in the regular
#' # simulate() call:
#' control$MCMC.samplesize <- 3
#'
#' # Capture intermediate networks; can also be left NULL for just the
#' # statistics:
#' control$MCMC.save_networks <- TRUE
#'
#' # Simulate starting from this state:
#' out <- ergm_MCMC_sample(sim.state$object, control, theta = -1, verbose=6)
#' names(out)
#' out$stats # Sampled statistics
#' str(out$networks, 1) # Updated ergm_state (one per thread)
#' # List (an element per thread) of lists of captured ergm_states,
#' # one for each sampled network:
#' str(out$sampnetworks, 2)
#' lapply(out$sampnetworks[[1]], as.network) # Converted to networks.
#'
#' # One more, picking up where the previous sampler left off, but see Note:
#' control$MCMC.samplesize <- 1
#' str(ergm_MCMC_sample(out$networks, control, theta = -1, verbose=6), 2)
#'
#' @export
ergm_MCMC_sample <- function(state, control, theta=NULL, 
                             verbose=FALSE,..., eta=ergm.eta(theta, (if(is.ergm_state(state))as.ergm_model(state)else as.ergm_model(state[[1]]))$etamap)){
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, verbose)
  
  if(is.ergm_state(state)) state <- list(state)
  state <- rep(state, length.out=nthreads(control))

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

  handle_statuses <- function(outl){
    statuses <- map_int(outl,"status")
    if(2L %in% statuses) stop("Sampling failed due to a Metropolis-Hastings proposal failing.") # MCMC_MH_FAILED: MH proposal failed somewhere. Throw an error.
    else if(1L %in% statuses) 1L # MCMC_TOO_MANY_EDGES, exceeding even control.parallel$MCMC.maxedges
    else 0L
  }

  sms <- vector("list", nthreads(control))
  nws <- if(NVL(control$MCMC.save_networks, FALSE)) vector("list", nthreads(control))

  if(!is.null(control.parallel$MCMC.effectiveSize)){
    #################################
    ######### Adaptive MCMC #########
    #################################
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
        if(is.na(burnin.pval) | burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
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
      if(status <- handle_statuses(outl)) return(list(status=status)) # Stop if something went wrong.

      sms <- mapply(rbind, sms, map(outl, "s"), SIMPLIFY=FALSE)
      if(!is.null(nws)) nws <- mapply(c, nws, map(outl, "saved"), SIMPLIFY=FALSE)
      state <- map(outl, "state")
      
      while(nrow(sms[[1]])-best.burnin$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)){
          sms[[i]] <- sms[[i]][seq_len(floor(nrow(sms[[i]])/2))*2+nrow(sms[[i]])%%2,,drop=FALSE]
          if(!is.null(nws)) nws[[i]] <- nws[[i]][seq_len(floor(length(nws[[i]])/2))*2+length(nws[[i]])%%2]
        }
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
      if(is.na(best.burnin$burnin)){
        if(verbose>1) message("Can not compute a valid burn-in. Setting burn-in to",interval,".")
        best.burnin$burnin <- interval
      }
      burnin.pval <- best.burnin$pval
      if(is.na(burnin.pval) | burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
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
      if(!is.null(nws)) nws[[i]] <- nws[[i]][seq_len(floor(length(nws[[i]])/2))*2+length(nws[[i]])%%2]
      outl[[i]]$final.interval <- interval
    }
  }else{
    #################################
    ########## Static MCMC ##########
    #################################
    outl <- doruns()
    if(status <- handle_statuses(outl)) return(list(status=status)) # Stop if something went wrong.
    sms <- map(outl, "s") %>% map(coda::mcmc, control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
    if(!is.null(nws)) nws <- map(outl, "saved")
    
    if(control.parallel$MCMC.runtime.traceplot){
      lapply(sms, function(sm) NVL3(theta, ergm.estfun(sm, ., as.ergm_model(state[[1]])), sm[,!as.ergm_model(state[[1]])$etamap$offsetmap,drop=FALSE])) %>% lapply(mcmc, start=control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval) %>% as.mcmc.list() %>% window(., thin=thin(.)*max(1,floor(niter(.)/1000))) %>% plot(ask=FALSE,smooth=TRUE,density=FALSE)
    }
  }

  #
  #   Process the results
  #
  statsmatrices <- vector("list", nthreads(control))
  newnetworks <- vector("list", nthreads(control))
  sampnetworks <- if(!is.null(nws)) vector("list", nthreads(control))

  final.interval <- c()
  for(i in (1:nthreads(control))){
    z <- outl[[i]]
    statsmatrices[[i]] <- sms[[i]]
    newnetworks[[i]] <- update(state0[[i]], state=z$state)
    if(!is.null(nws)) sampnetworks[[i]] <- lapply(nws[[i]], function(state) update(state0[[i]], state=state))
    final.interval <- c(final.interval, z$final.interval)
  }
  
  stats <- as.mcmc.list(statsmatrices)
  if(verbose){message("Sample size = ",niter(stats)*nchain(stats)," by ",
                  niter(stats),".")}
  
  list(stats = stats, networks=newnetworks, sampnetworks=sampnetworks, status=0, final.interval=final.interval)
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

  MCMC.maxedges <- NVL(control$MCMC.maxedges, Inf)
  if(NVL(control$MCMC.save_networks, FALSE)) MCMC.maxedges <- -MCMC.maxedges

  z <-
    if(!is.valued(state))
      .Call("MCMC_wrapper",
            state,
            # MCMC settings
            as.double(deInf(eta)),
            as.integer(samplesize),
            as.integer(burnin), 
            as.integer(interval),
            as.integer(deInf(MCMC.maxedges, "maxint")),
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
            as.integer(deInf(MCMC.maxedges, "maxint")),
            as.integer(verbose),
            PACKAGE="ergm")

  if(z$status) return(z) # If there is an error.
  z$s <- matrix(z$s, ncol=nparam(state,canonical=TRUE), byrow = TRUE)
  colnames(z$s) <- param_names(state, canonical=TRUE)
  z$state <- ergm_state_receive(z$state)
  z$saved <- EVL(lapply(z$saved, ergm_state_receive))

  z
}


.find_OK_burnin <- function(x, ...){
  n <- nrow(x[[1]])
  ssr <- function(b, s){
    b <- round(b)
    a <- lm(s ~ c(seq_len(b) - 1, rep(b, n - b)))
    sum(sigma(a)^2)
  }
  geweke <- function(b){
    if(b>0) x <- window(x, start=start(x) + b * thin(x))
    p.val <- suppressWarnings(geweke.diag.mv(x, ...)$p.value)
    if(is.na(p.val)) 0 else p.val
  }

  FAIL <- list(burnin=round(n/2), pval=0)
  xs <- x %>% map(scale) %>% map(~.[,attr(.,"scaled:scale")>0,drop=FALSE]) %>% discard(~ncol(.)==0)
  if(length(xs)==0) return(FAIL)

  best <- sapply(xs, function(x) optimize(ssr, c(0, n/2), s=x, tol=1)$minimum)
  if(all(is.na(best) | is.infinite(best))) return(FAIL)

  best <- max(best, na.rm=TRUE)
  list(burnin=round(best), pval=geweke(round(best)))
}
