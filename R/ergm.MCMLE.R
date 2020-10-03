#  File R/ergm.MCMLE.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
############################################################################
# The <ergm.MCMLE> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is the default and uses
# optimization of an MCMC estimate of the log-likelihood.  (The other
# MLE styles are found in functions <ergm.robmon>, <ergm.stocapprox>, and
# <ergm.stepping> 
#
#
# --PARAMETERS--
#   init         : the initial theta values
#   nw             : the network 
#   model          : the model, as returned by <ergm_model>
#   initialfit     : an ergm object, as the initial fit, possibly returned
#                    by <ergm.initialfit>
#   control     : a list of parameters for controlling the MCMC sampling;
#                    recognized components include
#       samplesize : the number of MCMC sampled networks
#       maxit      : the maximum number of iterations to use
#       parallel   : the number of threads in which to run the sampling
#       packagenames: names of packages; this is only relevant if "ergm" is given
#       interval    : the number of proposals to ignore between sampled networks
#       burnin      : the number of proposals to initially ignore for the burn-in
#                     period
#
#       epsilon    : ??, this is essentially unused, except to print it if
#                    'verbose'=T and to pass it along to <ergm.estimate>,
#                    which ignores it;   
#   proposal     : an proposal object for 'nw', as returned by
#                    <proposal>
#   proposal.obs : an proposal object for the observed network of'nw',
#                    as returned by <proposal>
#   verbose        : whether the MCMC sampling should be verbose (T or F);
#                    default=FALSE
#   sequential     : whether to update the network returned in
#                    'v$newnetwork'; if the network has missing edges,
#                    this is ignored; default=control$MCMLE.sequential
#   estimate       : whether to optimize the init coefficients via
#                    <ergm.estimate>; default=TRUE
#   ...            : additional parameters that may be passed from within;
#                    all are ignored
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.MCMLE>=*);
#      note that if the model is degenerate, only 'coef' and 'sample' are
#      returned; if 'estimate'=FALSE, the MCMC and se variables will be
#      NA or NULL
#
#############################################################################

ergm.MCMLE <- function(init, nw, model,
                             initialfit, 
                             control, 
                             proposal, proposal.obs,
                             verbose=FALSE,
                             sequential=control$MCMLE.sequential,
                             estimate=TRUE,
                             response=NULL, ...) {
  message("Starting Monte Carlo maximum likelihood estimation (MCMLE):")
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  steplen.hist <- c()
  steplen <- control$MCMLE.steplength
  if(control$MCMLE.steplength=="adaptive") steplen <- 1

  control$MCMC.effectiveSize <- control$MCMLE.effectiveSize
  control$MCMC.base.samplesize <- control$MCMC.samplesize
  control$obs.MCMC.base.samplesize <- control$obs.MCMC.samplesize

  control0 <- control

  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, max(verbose-1,0))
  
  # Store information about original network, which will be returned at end
  nw.orig <- nw

  # Impute missing dyads.
  #
  # Note: We do not need to update nw.stats, because if we are in a
  # situation where we are imputing dyads, the optimization is in the
  # observational mode, and since both the constrained and the
  # unconstrained samplers start from the same place, the initial
  # statshifts will be 0. target.stats and missing dyads are mutually
  # exclusive, so model$target.stats will be set equal to
  # model$nw.stats, causing this to happen.
  nw <- single.impute.dyads(nw, response=response, constraints=proposal$arguments$constraints, constraints.obs=proposal.obs$arguments$constraints, min_informative = control$obs.MCMC.impute.min_informative, default_density = control$obs.MCMC.impute.default_density, output="pending", verbose=verbose)

  ec <- if(is(nw, "network")) network.edgecount(nw, FALSE)
        else nrow(as.edgelist(nw))

  if(control$MCMLE.density.guard>1){
    # Calculate the density guard threshold.
    control$MCMC.max.maxedges <- round(min(control$MCMC.max.maxedges,
                                           max(control$MCMLE.density.guard*ec,
                                               control$MCMLE.density.guard.min)))
    control$MCMC.init.maxedges <- round(min(control$MCMC.max.maxedges, control$MCMC.init.maxedges))
    if(verbose) message("Density guard set to ",control$MCMC.max.maxedges," from an initial count of ",ec," edges.")
  }  

  nws <- rep(list(nw),nthreads(control)) # nws is now a list of networks.

  # statshift is the difference between the target.stats (if
  # specified) and the statistics of the networks in the LHS of the
  # formula or produced by SAN. If target.stats is not speficied
  # explicitly, they are computed from this network, so
  # statshift==0. To make target.stats play nicely with offsets, we
  # set statshifts to 0 where target.stats is NA (due to offset).
  statshift <- model$nw.stats - model$target.stats
  statshift[is.na(statshift)] <- 0
  statshifts <- rep(list(statshift), nthreads(control)) # Each network needs its own statshift.

  # Is there observational structure?
  obs <- ! is.null(proposal.obs)
  
  # Initialize control.obs and other *.obs if there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$MCMC.base.samplesize <- control$obs.MCMC.base.samplesize
    control.obs$MCMC.samplesize <- control$obs.MCMC.samplesize
    control.obs$MCMC.interval <- control$obs.MCMC.interval
    control.obs$MCMC.burnin <- control$obs.MCMC.burnin
    control.obs$MCMC.burnin.min <- control$obs.MCMC.burnin.min
    control0.obs <- control.obs

    nws.obs <- lapply(nws, identity)
    statshifts.obs <- statshifts
  }
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  calc.MCSE <- FALSE
  last.adequate <- FALSE

  # STATE_VARIABLES = variables collectively containing the state of the optimizer that would allow it to resume, excluding control lists.
  # CONTROL_VARIABLES = control lists
  # INTERMEDIATE_VARIABLES = variables of interest in debugging and diagnostics.
  #
  # Both lists need to be kept up to date with the implementation.
  STATE_VARIABLES <- c("nws", "nws.obs", "mcmc.init", "statshifts", "statshifts.obs", "calc.MCSE", "last.adequate", "coef.hist", "stats.hist", "stats.obs.hist", "steplen.hist", "steplen")
  CONTROL_VARIABLES <- c("control", "control.obs", "control0", "control0.obs")
  INTERMEDIATE_VARIABLES <- c("nws", "nws.obs", "statsmatrices", "statsmatrices.obs", "statshifts", "statshifts.obs", "coef.hist", "stats.hist", "stats.obs.hist", "steplen.hist")

  if(!is.null(control$resume)){
    message("Resuming from state saved in ", sQuote(control$resume),".")
    state <- new.env()
    load(control$resume,envir=state)

    # Merge control lists intelligently:
    .merge_controls <- function(saved.ctrl, saved.ctrl0, ctrl0){
      ctrl <- saved.ctrl
      for(name in union(names(ctrl), names(ctrl0))){
        if(!identical(ctrl[[name]],ctrl0[[name]])){ # Settings differ.
          if(!identical(ctrl0[[name]], saved.ctrl0[[name]])){
            if(verbose) message("Passed-in control setting ", sQuote(name), " changed from original run: overriding saved state.")
            ctrl[[name]] <- ctrl0[[name]]
          }else if(verbose) message("Passed-in control setting ", sQuote(name), " unchanged from original run: using saved state.")
        }
      }
      ctrl
    }

    control <- .merge_controls(state$control, state$control0, control0)
    if(obs) control.obs <- .merge_controls(state$control.obs, state$control0.obs, control0.obs)

    # Copy the rest
    for(name in intersect(ls(state), STATE_VARIABLES)) assign(name, state[[name]])

    # Clean up
    rm(state)
  }

  
  for(iteration in 1:control$MCMLE.maxit){
    if(verbose){
      message("\nIteration ",iteration," of at most ", control$MCMLE.maxit,
          " with parameter:")
      message_print(mcmc.init)
    }else{
      message("Iteration ",iteration," of at most ", control$MCMLE.maxit,":")
    }

    if(!is.null(control$checkpoint)){
      save(list=intersect(ls(), c(STATE_VARIABLES, CONTROL_VARIABLES)), file=sprintf(control$checkpoint, iteration))
    }

    # Obtain MCMC sample
    if(verbose) message("Starting unconstrained MCMC...")
    z <- ergm_MCMC_sample(nws, model, proposal, control, theta=mcmc.init, response=response, update.nws=FALSE, verbose=max(verbose-1,0))
        
    if(z$status==1) stop("Number of edges in a simulated network exceeds that in the observed by a factor of more than ",floor(control$MCMLE.density.guard),". This is a strong indicator of model degeneracy or a very poor starting parameter configuration. If you are reasonably certain that neither of these is the case, increase the MCMLE.density.guard control.ergm() parameter.")
        
    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrices <- mapply(sweep, z$stats, statshifts, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
    for(i in seq_along(statsmatrices)) colnames(statsmatrices[[i]]) <- param_names(model,canonical=TRUE)
    nws.returned <- z$networks
    statsmatrix <- do.call(rbind,statsmatrices)
    
    if(verbose){
      message("Back from unconstrained MCMC.")
      if(verbose>1){
        message("Average statistics:")
        message_print(colMeans(statsmatrix))
      }
    }
    
    ##  Does the same, if observation process:
    if(obs){
      if(verbose) message("Starting constrained MCMC...")
      z.obs <- ergm_MCMC_sample(nws.obs, model, proposal.obs, control.obs, theta=mcmc.init, response=response, update.nws=FALSE, verbose=max(verbose-1,0))
      
      if(z.obs$status==1) stop("Number of edges in the simulated network exceeds that observed by a large factor (",control$MCMC.max.maxedges,"). This is a strong indication of model degeneracy. If you are reasonably certain that this is not the case, increase the MCMLE.density.guard control.ergm() parameter.")
      
      statsmatrices.obs <- mapply(sweep, z.obs$stats, statshifts.obs, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
      for(i in seq_along(statsmatrices.obs)) colnames(statsmatrices.obs[[i]]) <- param_names(model,canonical=TRUE)
      nws.obs.returned <- z.obs$networks
      statsmatrix.obs <- do.call(rbind,statsmatrices.obs)
      
      if(verbose){
        message("Back from constrained MCMC.")
        if(verbose>1){
          message("Average statistics:")
          message_print(colMeans(statsmatrix.obs))
        }
      }
    }else{
      statsmatrices.obs <- statsmatrix.obs <- NULL
      z.obs <- NULL
    }
    
    if(sequential) {
      nws <- nws.returned
      statshifts <- lapply(statsmatrices, function(sm) sm[nrow(sm),])
      
      if(obs){
        nws.obs <- nws.obs.returned
        statshifts.obs <- lapply(statsmatrices.obs, function(sm) sm[nrow(sm),])
      }      
    }

    if(!is.null(control$MCMLE.save_intermediates)){
      save(list=intersect(ls(), INTERMEDIATE_VARIABLES), file=sprintf(control$MCMLE.save_intermediates, iteration))
    }

    # Compute the sample estimating functions and the convergence p-value. 
    esteq <- ergm.estfun(statsmatrix, mcmc.init, model)
    if(isTRUE(all.equal(apply(esteq,2,stats::sd), rep(0,ncol(esteq)), check.names=FALSE))&&!all(esteq==0))
      stop("Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.")

    check_nonidentifiability(esteq, NULL, model,
                             tol = control$MCMLE.nonident.tol, type="statistics",
                             action = control$MCMLE.nonident)

    esteq.obs <- if(obs) ergm.estfun(statsmatrix.obs, mcmc.init, model) else NULL

    # Update the interval to be used.
    if(!is.null(control$MCMC.effectiveSize)){
      control$MCMC.interval <- round(max(z$final.interval/control$MCMLE.effectiveSize.interval_drop,1))
      if(verbose) message("New interval = ",control$MCMC.interval,".")
      if(obs){
        control.obs$MCMC.interval <- round(max(z.obs$final.interval/control$MCMLE.effectiveSize.interval_drop,1))
        if(verbose) message("New constrained interval = ",control.obs$MCMC.interval,".")
      }
    }
        
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..

    # These are only nontrivial when the model is curved or when there are missing data.
    if(verbose){
      message("Average estimating function values:")
      message_print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else -colMeans(esteq))
    }

    if(!estimate){
      if(verbose){message("Skipping optimization routines...")}
      nws.returned <- lapply(nws.returned, as.network, response=response)
      l <- list(coef=mcmc.init, mc.se=rep(NA,length=length(mcmc.init)),
                sample=statsmatrix, sample.obs=statsmatrix.obs,
                iterations=1, MCMCtheta=mcmc.init,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(mcmc.init)), #acf=NULL,
                samplesize=control$MCMC.samplesize, failure=TRUE,
                newnetwork = nws.returned[[1]],
                newnetworks = nws.returned)
      return(structure (l, class="ergm"))
    } 

    statsmatrix.0 <- statsmatrix
    statsmatrix.0.obs <- statsmatrix.obs
    if(control$MCMLE.steplength=="adaptive"){
      if(verbose){message("Starting adaptive MCMLE Optimization...")}
      adaptive.steplength <- 2
      statsmean <- apply(statsmatrix.0,2,base::mean)
      v <- list(loglikelihood=control$MCMLE.adaptive.trustregion*2)
      while(v$loglikelihood > control$MCMLE.adaptive.trustregion){
        adaptive.steplength <- adaptive.steplength / 2
        if(!is.null(statsmatrix.0.obs)){
          statsmatrix.obs <- t(adaptive.steplength*t(statsmatrix.0.obs) + (1-adaptive.steplength)*statsmean) # I.e., shrink each point of statsmatrix.obs towards the centroid of statsmatrix.
        }else{
          statsmatrix <- sweep(statsmatrix.0,2,(1-adaptive.steplength)*statsmean,"-")
        }
        if(verbose){message("Optimizing with step length ",adaptive.steplength,".")}
        #
        #   If not the last iteration do not compute all the extraneous
        #   statistics that are not needed until output
        #
        v<-ergm.estimate(init=mcmc.init, model=model,
                         statsmatrix=statsmatrix, 
                         statsmatrix.obs=statsmatrix.obs, 
                         epsilon=control$epsilon,
                         nr.maxit=control$MCMLE.NR.maxit,
                         nr.reltol=control$MCMLE.NR.reltol,
                         calc.mcmc.se=control$MCMLE.termination == "precision" || (control$MCMC.addto.se && last.adequate) || iteration == control$MCMLE.maxit, 
                         hessianflag=control$main.hessian,
                         trustregion=control$MCMLE.trustregion, method=control$MCMLE.method,
                         metric=control$MCMLE.metric,
                         dampening=control$MCMLE.dampening,
                         dampening.min.ess=control$MCMLE.dampening.min.ess,
                         dampening.level=control$MCMLE.dampening.level,
                         compress=control$MCMC.compress, verbose=verbose,
                         estimateonly=TRUE)
      }
      if(v$loglikelihood < control$MCMLE.trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        message("The log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),".")
        options(scipen=current.scipen)
      }else{
        message("The log-likelihood did not improve.")
      }
      steplen.hist <- c(steplen.hist, adaptive.steplength)
      steplen <- adaptive.steplength
    }else{
      if(verbose){message("Starting MCMLE Optimization...")}
      steplen <-
        if(!is.null(control$MCMLE.steplength.margin))
          .Hummel.steplength(
            if(control$MCMLE.steplength.esteq) esteq else statsmatrix.0[,!model$etamap$offsetmap,drop=FALSE], 
            if(control$MCMLE.steplength.esteq) esteq.obs else statsmatrix.0.obs[,!model$etamap$offsetmap,drop=FALSE],
            control$MCMLE.steplength.margin, control$MCMLE.steplength,steplength.prev=steplen,verbose=verbose,
            x2.num.max=control$MCMLE.steplength.miss.sample, steplength.maxit=control$MCMLE.steplength.maxit,
            parallel=control$MCMLE.steplength.parallel, control=control
          )
        else control$MCMLE.steplength
      
      if(steplen==control$MCMLE.steplength || is.null(control$MCMLE.steplength.margin) || iteration==control$MCMLE.maxit) calc.MCSE <- TRUE
      
      statsmean <- apply(statsmatrix.0,2,base::mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- t(steplen*t(statsmatrix.0.obs) + (1-steplen)*statsmean) # I.e., shrink each point of statsmatrix.obs towards the centroid of statsmatrix.
      }else{
        statsmatrix <- sweep(statsmatrix.0,2,(1-steplen)*statsmean,"-")
      }
      steplen.hist <- c(steplen.hist, steplen)
      
      message("Optimizing with step length ",steplen,".")
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(init=mcmc.init, model=model,
                       statsmatrix=statsmatrix, 
                       statsmatrix.obs=statsmatrix.obs, 
                       epsilon=control$epsilon,
                       nr.maxit=control$MCMLE.NR.maxit,
                       nr.reltol=control$MCMLE.NR.reltol,
                       calc.mcmc.se=control$MCMLE.termination == "precision" || (control$MCMC.addto.se && last.adequate) || iteration == control$MCMLE.maxit,
                       hessianflag=control$main.hessian,
                       trustregion=control$MCMLE.trustregion, 
                       method=control$MCMLE.method,
                       dampening=control$MCMLE.dampening,
                       dampening.min.ess=control$MCMLE.dampening.min.ess,
                       dampening.level=control$MCMLE.dampening.level,
                       metric=control$MCMLE.metric,
                       compress=control$MCMC.compress, verbose=verbose,
                       estimateonly=!calc.MCSE)
      if(v$loglikelihood < control$MCMLE.trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        message("The log-likelihood improved by ",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),".")
        options(scipen=current.scipen)
      }else{
        message("The log-likelihood did not improve.")
      }
    }
          
    mcmc.init <- v$coef
    coef.hist <- rbind(coef.hist, mcmc.init)
    stats.obs.hist <- NVL3(statsmatrix.obs, rbind(stats.obs.hist, apply(.[], 2, base::mean)))
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, base::mean))

    # This allows premature termination.
    
    if(steplen<control$MCMLE.steplength){ # If step length is less than its maximum, don't bother with precision stuff.
      last.adequate <- FALSE
      control$MCMC.samplesize <- control$MCMC.base.samplesize
      if(obs) control.obs$MCMC.samplesize <- control.obs$MCMC.base.samplesize
      
    } else {
    
    if(control$MCMLE.termination == "precision"){
      prec.loss <- (sqrt(diag(v$mc.cov+v$covar))-sqrt(diag(v$covar)))/sqrt(diag(v$mc.cov+v$covar))
      if(verbose){
        message("Standard Error:")
        message_print(sqrt(diag(v$covar)))
        message("MC SE:")
        message_print(sqrt(diag(v$mc.cov)))
        message("Linear scale precision loss due to MC estimation of the likelihood:")
        message_print(prec.loss)
      }
      if(sqrt(mean(prec.loss^2, na.rm=TRUE)) <= control$MCMLE.MCMC.precision){
        if(last.adequate){
          message("Precision adequate twice. Stopping.")
          break
        }else{
          message("Precision adequate. Performing one more iteration.")
          last.adequate <- TRUE
        }
      }else{
        last.adequate <- FALSE
        prec.scl <- max(sqrt(mean(prec.loss^2, na.rm=TRUE))/control$MCMLE.MCMC.precision, 1) # Never decrease it.
        
        if (!is.null(control$MCMC.effectiveSize)) { # ESS-based sampling
          control$MCMC.effectiveSize <- round(control$MCMC.effectiveSize * prec.scl)
          if(control$MCMC.effectiveSize/control$MCMC.samplesize>control$MCMLE.MCMC.max.ESS.frac) control$MCMC.samplesize <- control$MCMC.effectiveSize/control$MCMLE.MCMC.max.ESS.frac
          # control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
          message("Increasing target MCMC sample size to ", control$MCMC.samplesize, ", ESS to",control$MCMC.effectiveSize,".")
        } else { # Fixed-interval sampling
          control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
          control$MCMC.burnin <- round(control$MCMC.burnin * prec.scl)
          message("Increasing MCMC sample size to ", control$MCMC.samplesize, ", burn-in to",control$MCMC.burnin,".")
        }

        if(obs){
          if (!is.null(control$MCMC.effectiveSize)) { # ESS-based sampling
            control$MCMC.effectiveSize <- round(control$MCMC.effectiveSize * prec.scl)
            if(control$MCMC.effectiveSize/control.obs$MCMC.samplesize>control$MCMLE.MCMC.max.ESS.frac) control.obs$MCMC.samplesize <- control$MCMC.effectiveSize/control$MCMLE.MCMC.max.ESS.frac
            # control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
            message("Increasing target constrained MCMC sample size to ", control.obs$MCMC.samplesize, ", ESS to",control$MCMC.effectiveSize,".")
          } else { # Fixed-interval sampling
            control.obs$MCMC.samplesize <- round(control.obs$MCMC.samplesize * prec.scl)
            control.obs$MCMC.burnin <- round(control.obs$MCMC.burnin * prec.scl)
            message("Increasing constrained MCMC sample size to ", control.obs$MCMC.samplesize, ", burn-in to",control.obs$MCMC.burnin,".")
          }
        }
      }
    }else if(control$MCMLE.termination=='Hotelling'){
      conv.pval <- ERRVL(try(approx.hotelling.diff.test(esteq, esteq.obs)$p.value), NA)
      message("Nonconvergence test p-value:",conv.pval,"")
      if(!is.na(conv.pval) && conv.pval>=control$MCMLE.conv.min.pval){
        message("No nonconvergence detected. Stopping.")
        break
      }      
    }else if(control$MCMLE.termination=='Hummel'){
      if(last.adequate){
        message("Step length converged twice. Stopping.")
        break
      }else{
        message("Step length converged once. Increasing MCMC sample size.")
        last.adequate <- TRUE
        control$MCMC.samplesize <- control$MCMC.base.samplesize * control$MCMLE.last.boost
        if(obs) control.obs$MCMC.samplesize <- control.obs$MCMC.base.samplesize * control$MCMLE.last.boost
      }
    }
    
    }

    #' @importFrom utils tail
    # stop if MCMLE is stuck (steplen stuck near 0)
    if ((length(steplen.hist) > 2) && sum(tail(steplen.hist,2)) < 2*control$MCMLE.steplength.min) {
      stop("MCMLE estimation stuck. There may be excessive correlation between model terms, suggesting a poor model for the observed data. If target.stats are specified, try increasing SAN parameters.")
    }    
    #Otherwise, don't stop before iterations are exhausted.
    if (iteration == control$MCMLE.maxit) {
      message("MCMLE estimation did not converge after ", control$MCMLE.maxit, " iterations. The estimated coefficients may not be accurate. Estimation may be resumed by passing the coefficients as initial values; see 'init' under ?control.ergm for details.")
    }
  } # end of main loop

  message("Finished MCMLE.")

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, control)

  nws.returned <- lapply(nws.returned, as.network, response=response)
  v$network <- nw.orig
  v$newnetworks <- nws.returned
  v$newnetwork <- nws.returned[[1]]
  v$coef.init <- init
  #v$initialfit <- initialfit
  v$est.cov <- v$mc.cov
  v$mc.cov <- NULL

  v$coef.hist <- coef.hist
  v$stats.hist <- stats.hist
  v$stats.obs.hist <- stats.obs.hist
  v$steplen.hist <- steplen.hist
  
  v$iterations <- iteration
  v$control <- control
  
  v$etamap <- model$etamap
  v
}

