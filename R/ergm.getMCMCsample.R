#  File R/ergm.getMCMCsample.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
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
#' \item{final.effectiveSize}{adaptively determined target ESS (non-trivial if `control$MCMC.effectiveSize` is specified via a matrix).}
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

  utils::flush.console()

  force(eta)
  state0 <- state

  send_model_proposal <- function(){
    if(verbose>1) message("Populating the state cache on worker nodes.")
    call_state_cache <- function(...) ergm::ergm_state_cache(...)
    environment(call_state_cache) <- globalenv()

    for(item in c("model", "proposal")){
      have_item <- unlist(clusterMap(ergm.getCluster(control), call_state_cache,
                                     list("check"), map(state0, c("uids", item))))
      if(verbose>1 && all(have_item)) message("State cache for ", item, " already populated.")
      if(!all(have_item))
        clusterMap(ergm.getCluster(control), call_state_cache,
                   ifelse(have_item, "pass", "insert"), map(state0, c("uids", item)), ifelse(have_item, list(NULL), map(state0, item)))
    }
  }

  if(!is.null(ergm.getCluster(control))){ # Populate cache and
    send_model_proposal()
    state <- lapply(state, ergm_state_receive) # Don't carry around nw0, model, or proposal.
    # Don't carry the environment with the function.
    call_MCMC_worker <- function(...) ergm::ergm_MCMC_slave(...)
    environment(call_MCMC_worker) <- globalenv()
  }else state <- lapply(state, ergm_state_send) # Don't carry around nw0.
  
  #' @importFrom parallel clusterMap
  doruns <- function(burnin=NULL, samplesize=NULL, interval=NULL){
    if(!is.null(ergm.getCluster(control)))
      persistEvalQ({
        clusterMap(ergm.getCluster(control), call_MCMC_worker,
                   state=state, MoreArgs=list(eta=eta,control=control.parallel,verbose=verbose,...,burnin=burnin,samplesize=samplesize,interval=interval))},
        retries = getOption("ergm.cluster.retries"),
        beforeRetry = {
          ergm.restartCluster(control,verbose)
          send_model_proposal()
        })
    else{
      out <- list(ergm_MCMC_slave(state[[1]], burnin=burnin,samplesize=samplesize,interval=interval,eta=eta,control=control.parallel,verbose=verbose,...))
      # Note: the return value's state will be a ergm_state_receive.
      for(i in seq_along(out)) out[[i]]$state <- update(state[[i]], out[[i]]$state)
      out
    }
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
    if(verbose) message("Starting adaptive MCMC with target ESS ", format(control.parallel$MCMC.effectiveSize), "...")

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
          if(verbose>1) message("Predicted additional sample size: ", format(pred.ss), " dampened to ", format(damp.ss), ", so running ", samplesize, " steps forward.")
        }
      }
        
      outl<-doruns(burnin = interval, # I.e., skip that much before the first draw.
                   samplesize = samplesize,
                   interval = interval)
      if(status <- handle_statuses(outl)) return(list(status=status)) # Stop if something went wrong.

      sms <- Map(rbind, sms, map(outl, "s"))
      if(!is.null(nws)) nws <- Map(c, nws, map(outl, "saved"))
      state <- map(outl, "state")
      
      while(nrow(sms[[1]])-best.burnin$burnin>=(control.parallel$MCMC.samplesize)*2){
        for(i in seq_along(outl)){
          sms[[i]] <- sms[[i]][seq_len(floor(nrow(sms[[i]])/2))*2+nrow(sms[[i]])%%2,,drop=FALSE]
          if(!is.null(nws)) nws[[i]] <- nws[[i]][seq_len(floor(length(nws[[i]])/2))*2+length(nws[[i]])%%2]
        }
        interval <- interval*2
        if(verbose) message("Increasing thinning to ",interval,".")
      }
      
      esteq <- lapply(sms, function(sm) NVL3(theta, ergm.estfun(sm, ., as.ergm_model(state0[[1]])), sm[,!as.ergm_model(state0[[1]])$etamap$offsetmap,drop=FALSE])) %>%
        lapply.mcmc.list(mcmc, start=1, thin=interval) %>% lapply.mcmc.list(`-`)

      if(control.parallel$MCMC.runtime.traceplot){
        plot(window(esteq, thin=thin(esteq)*max(1,floor(niter(esteq)/1000)))
            ,ask=FALSE,smooth=TRUE,density=FALSE)
      }

      best.burnin <- .find_OK_burnin(esteq, control)
      burnin.pval <- best.burnin$pval
      postburnin.mcmc <- window(esteq, start=start(esteq)+best.burnin$burnin*thin(esteq))

      if(is.na(burnin.pval) || burnin.pval <= control$MCMC.effectiveSize.burnin.pval){
        if(is.const.sample(postburnin.mcmc)){
          message("Post-burnin sample is constant; returning.")
          control.parallel$MCMC.effectiveSize <- eS <- 1
          break
        }
        if(verbose>1) message("Selected burn-in ", format(start(esteq)+best.burnin$burnin*thin(esteq), digits=2, scientific=TRUE), " (",round(best.burnin$burnin/niter(esteq)*100,2),"%) p-value = ", format(burnin.pval), " is below the threshold of ",control$MCMC.effectiveSize.burnin.pval,".")
        next
      }

      if(is.matrix(control$MCMC.effectiveSize)){
        ## Target ESS determined by covariance matrix.
        postburnin.var <- var.mcmc.list(postburnin.mcmc)
        control.parallel$MCMC.effectiveSize <- mean(diag(control$MCMC.effectiveSize%*%postburnin.var))
        if(verbose) message("Variance-based adaptive MCMC set target ESS to ", format(control.parallel$MCMC.effectiveSize), ".")
      }
      
      eS <- niter(postburnin.mcmc)*nchain(postburnin.mcmc)/attr(spectrum0.mvar(postburnin.mcmc, order.max=control$MCMC.effectiveSize.order.max),"infl")
      
      if(verbose) message("ESS of ", format(eS)," attained with burn-in of ", round(best.burnin$burnin/niter(esteq)*100,2),"%; convergence p-value = ", format(burnin.pval), ".")

      if(eS>=control.parallel$MCMC.effectiveSize){
        if(burnin.pval > control$MCMC.effectiveSize.burnin.pval){
          if(verbose) message("Target ESS achieved and is trustworthy. Returning.")
          break
        }else{
          if(verbose>1) message("ESS and burn-in estimates are not trustworthy.")
        }
      }
    }

    if(control.parallel$MCMC.effectiveSize > 1 && eS < control.parallel$MCMC.effectiveSize)
      warning("Unable to reach target ESS in iterations alotted.")

    for(i in seq_along(outl)){
      if(best.burnin$burnin) sms[[i]] <- sms[[i]][-seq_len(best.burnin$burnin),,drop=FALSE]
      sms[[i]] <- coda::mcmc(sms[[i]], (best.burnin$burnin+1)*interval, thin=interval)
      if(!is.null(nws)) nws[[i]] <- nws[[i]][seq_len(floor(length(nws[[i]])/2))*2+length(nws[[i]])%%2]
    }

    final.interval <- interval
    final.effectiveSize <- control.parallel$MCMC.effectiveSize
  }else{
    #################################
    ########## Static MCMC ##########
    #################################
    outl <- doruns()
    if(status <- handle_statuses(outl)) return(list(status=status)) # Stop if something went wrong.
    sms <- map(outl, "s") %>% map(coda::mcmc, control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval)
    if(!is.null(nws)) nws <- map(outl, "saved")
    
    if(control.parallel$MCMC.runtime.traceplot){
      lapply(sms, function(sm) NVL3(theta, ergm.estfun(sm, ., as.ergm_model(state0[[1]])), sm[,!as.ergm_model(state0[[1]])$etamap$offsetmap,drop=FALSE])) %>% lapply(mcmc, start=control.parallel$MCMC.burnin+1, thin=control.parallel$MCMC.interval) %>% as.mcmc.list() %>% window(., thin=thin(.)*max(1,floor(niter(.)/1000))) %>% plot(ask=FALSE,smooth=TRUE,density=FALSE)
    }

    final.effectiveSize <- final.interval <- NULL
  }

  #
  #   Process the results
  #
  statsmatrices <- vector("list", nthreads(control))
  newnetworks <- vector("list", nthreads(control))
  sampnetworks <- if(!is.null(nws)) vector("list", nthreads(control))

  for(i in (1:nthreads(control))){
    z <- outl[[i]]
    statsmatrices[[i]] <- sms[[i]]
    newnetworks[[i]] <- update(state0[[i]], state=z$state)
    if(!is.null(nws)) sampnetworks[[i]] <- lapply(nws[[i]], function(state) update(state0[[i]], state=state))
  }
  
  stats <- as.mcmc.list(statsmatrices)
  if(verbose){message("Sample size = ",niter(stats)*nchain(stats)," by ",
                  niter(stats),".")}
  
  list(stats = stats, networks=newnetworks, sampnetworks=sampnetworks, status=0, final.interval=final.interval, final.effectiveSize=final.effectiveSize)
}

#' @rdname ergm_MCMC_sample
#' @description The \code{ergm_MCMC_slave} function calls the actual C
#'   routine and does minimal preprocessing.
#'
#' @param burnin,samplesize,interval MCMC paramters that can be used
#'   to temporarily override those in the `control` list.
#' @return \code{ergm_MCMC_slave} returns the MCMC sample as a list of
#'   the following: \item{s}{the matrix of statistics.}
#'   \item{state}{an [`ergm_state`] object for the new network.}
#'   \item{status}{success or failure code: `0` is success, `1` for
#'   too many edges, and `2` for a Metropolis-Hastings proposal failing,
#'   `-1` for [`ergm_model`] or [`ergm_proposal`] not passed and
#'   missing from the cache.}
#' @useDynLib ergm
#' @export
ergm_MCMC_slave <- function(state, eta,control,verbose,..., burnin=NULL, samplesize=NULL, interval=NULL){
  on.exit(ergm_Cstate_clear())
  state <- ergm_state_send(state)
  if(is.null(state$model) || is.null(state$proposal)) return(list(status=-1L))

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


.find_OK_burnin <- function(x, control){
  if(niter(x) < control$MCMC.effectiveSize.burnin.nmin) warning("The per-thread sample size for estimating burn-in is very small. This should probably be fixed in the calling function.")

  if((pc <- control$MCMC.effectiveSize.burnin.PC)>0 && pc<nvar(x)){
    v <- svd(scale(as.matrix(x)), nu=0, nv=pc)$v
    x <- lapply.mcmc.list(x, `%*%`, v)
  }

  nit <- niter(x)
  p <- nvar(x)
  bscl <- control$MCMC.effectiveSize.burnin.scl

  ssr <- function(decay, y, results=FALSE){
    # decay is basically the number of steps corresponding to halving of
    # the difference in the expected value of the variable at the
    # current MCMC draw from the ultimate expected value.

    x <- rep(2^(-seq_len(nit)/decay), length.out = NROW(y))
    a <- try(lm(y ~ x, x=results, y=results))
    if(results) structure(a, decay=decay) else sum(sigma(a)^2)
  }

  fit_decay <- function(y, interval){
    if(is.list(y)) y <- do.call(rbind, y)

    y <- y %>% scale %>% `[`(,attr(.,"scaled:scale")>0, drop=FALSE)

    if(ncol(y) == 0) return(NULL)

    decay <- optimize(ssr, interval, y=y)$minimum

    ssr(decay, y, results=TRUE)
  }

  geweke <- function(b){
    if(b>0) x <- window(x, start=start(x) + b * thin(x))
    if(niter(x)*nchain(x) > (nmax <- max(control$MCMC.effectiveSize.burnin.nmax, nvar(x)*2*4)))
    x <- lapply.mcmc.list(x, `[`, round(seq(from=1,to=niter(x),length.out=round(nmax/nchain(x)))), , drop=FALSE)

    p.val <- suppressWarnings(geweke.diag.mv(x, order.max=control$MCMC.effectiveSize.order.max)$p.value)
    if(is.na(p.val)) 0 else p.val
  }

  best_burnin <- function(coef, decay, s) -decay * log2(s/bscl/abs(coef))
  best_burnin.lm <- function(fit) best_burnin(if(is.matrix(coef(fit))) coef(fit)[2,] else coef(fit)[2], attr(fit,"decay"), sigma(fit))

  FAIL <- list(burnin=round(nit*control$MCMC.effectiveSize.burnin.max), pval=0)

  bestfits <- lapply(x, function(chain) lapply(seq_len(p), function(i) fit_decay(chain[,i], c(0,nit*log2(bscl)*8)))) %>% unlist(recursive=FALSE) %>% compact

  best <- ifelse(sapply(bestfits, function(fit) sd(resid(fit))/sd(fit$y)<1-1/bscl*2),
                 sapply(bestfits, best_burnin.lm),
                 round(nit*control$MCMC.effectiveSize.burnin.min))

  best <- sqrt(mean(pmin(pmax(best,0),nit*control$MCMC.effectiveSize.burnin.max)^2, na.rm=TRUE))

  if(all(is.na(best) | is.infinite(best))) return(FAIL)

  list(burnin=round(best), pval=geweke(round(best)))
}
