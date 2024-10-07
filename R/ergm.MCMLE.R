#  File R/ergm.MCMLE.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
############################################################################
# The <ergm.MCMLE> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is the default and uses
# optimization of an MCMC estimate of the log-likelihood.
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

ergm.MCMLE <- function(init, s, s.obs,
                             control, 
                             verbose=FALSE,
                             estimate=TRUE, ...) {
  message("Starting Monte Carlo maximum likelihood estimation (MCMLE):")
  # Is there observational structure?
  obs <- ! is.null(s.obs)
  model <- s$model
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, nparam(model, canonical=TRUE))
  stats.obs.hist <- matrix(NA, 0, nparam(model, canonical=TRUE))
  steplen.hist <- c()
  steplen <- control$MCMLE.steplength

  if(is.null(control$MCMLE.samplesize)) control$MCMLE.samplesize <- max(control$MCMLE.samplesize.min,control$MCMLE.samplesize.per_theta*nparam(model,canonical=FALSE, offset=FALSE))
  if(obs && is.null(control$obs.MCMLE.samplesize)) control$obs.MCMLE.samplesize <- max(control$obs.MCMLE.samplesize.min,control$obs.MCMLE.samplesize.per_theta*nparam(model,canonical=FALSE, offset=FALSE))

  control <- remap_algorithm_MCMC_controls(control, "MCMLE")
   
  control$MCMC.base.effectiveSize <- control$MCMC.effectiveSize
  control$obs.MCMC.base.effectiveSize <- control$obs.MCMC.effectiveSize
  
  control$MCMC.base.samplesize <- control$MCMC.samplesize
  control$obs.MCMC.base.samplesize <- control$obs.MCMC.samplesize

  adapt <- !is.null(control$MCMC.effectiveSize)
  adapt.obs.ESS <- obs && adapt && !is.null(control$obs.MCMC.effectiveSize)
  adapt.obs.var <- obs && adapt && !adapt.obs.ESS

  control0 <- control

  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, max(verbose-1,0))

  if(control$MCMLE.density.guard>1){
    # Calculate the density guard threshold.
    ec <- network.edgecount(s)
    control$MCMC.maxedges <- round(min(control$MCMC.maxedges,
                                       max(control$MCMLE.density.guard*ec,
                                           control$MCMLE.density.guard.min)))
    if(verbose) message("Density guard set to ",control$MCMC.maxedges," from an initial count of ",ec," edges.")
  }  

  s <- rep(list(s),nthreads(control)) # s is now a list of states.
  
  # Initialize control.obs and other *.obs if there is observation structure
  if(obs){
    control.obs <- control
    for(name in OBS_MCMC_CONTROLS) control.obs[[name]] <- control[[paste0("obs.", name)]]
    control0.obs <- control.obs

    s.obs <- rep(list(s.obs),nthreads(control))
  }

  ## A helper function to increase the MCMC sample size and target effective size by the specified factor.
  .boost_samplesize <- function(boost, base=FALSE){
    for(ctrl in c("control", if(obs && !adapt.obs.var) "control.obs")){
      control <- get(ctrl, parent.frame())
      sampsize.boost <-
        NVL2(control$MCMC.effectiveSize,
             boost^control$MCMLE.sampsize.boost.pow,
             boost)
    
      control$MCMC.samplesize <- round((if(base) control$MCMC.base.samplesize else control$MCMC.samplesize) * sampsize.boost)
      control$MCMC.effectiveSize <- NVL3((if(base) control$MCMC.base.effectiveSize else control$MCMC.effectiveSize), . * boost)
      
      control$MCMC.samplesize <- ceiling(max(control$MCMC.samplesize, control$MCMC.effectiveSize*control$MCMLE.min.depfac))
      
      assign(ctrl, control, parent.frame())
    }

    .set_obs_samplesize()
    NULL
  }

  ## A helper function to set the variance-based effective size and
  ## MCMC sample size for the constrained sample.  There is no great
  ## way to compute the corresponding raw sample size, so we'll just
  ## set it to the same as the unconstrained.
  .set_obs_samplesize <- function(){
    if(!adapt.obs.var) return()
    control.obs$MCMC.effectiveSize <<- sginv(esteq.var, tol=.Machine$double.eps^(3/4)) * control$MCMC.effectiveSize
    control.obs$MCMC.samplesize <<- max(control$obs.MCMLE.samplesize.min, ceiling(control$MCMC.samplesize * min(obs.ESS.adj * 1.2, 1))) # Fudge factor
    NULL
  }

  if(control$MCMLE.termination=='confidence'){
    estdiff.prev <- NULL
    d2.not.improved <- rep(FALSE, control$MCMLE.confidence.boost.lag)
  }

  
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  calc.MCSE <- FALSE
  last.adequate <- FALSE
  if(adapt.obs.var) obs.ESS.adj <- 1

  # ERGM_STATE_ELEMENTS = elements of the ergm_state objects (currently s and s.obs) that need to be saved.
  # STATE_VARIABLES = variables collectively containing the state of the optimizer that would allow it to resume, excluding control lists.
  # CONTROL_VARIABLES = control lists
  # INTERMEDIATE_VARIABLES = variables of interest in debugging and diagnostics.
  #
  # All lists need to be kept up to date with the implementation.
  ERGM_STATE_ELEMENTS <- c("el", "nw0", "stats", "ext.state", "ext.flag")
  STATE_VARIABLES <- c("mcmc.init", "calc.MCSE", "last.adequate", "coef.hist", "stats.hist", "stats.obs.hist", "steplen.hist", "steplen","setdiff.prev","d2.not.improved")
  CONTROL_VARIABLES <- c("control", "control.obs", "control0", "control0.obs")
  INTERMEDIATE_VARIABLES <- c("s", "s.obs", "statsmatrices", "statsmatrices.obs", "coef.hist", "stats.hist", "stats.obs.hist", "steplen.hist")

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

    # TODO: Implement a version with proper encapsulation.
    for(i in seq_along(s)) for(name in ERGM_STATE_ELEMENTS) s[[i]][[name]] <- state$s.reduced[[i]][[name]]
    if(obs) for(i in seq_along(s.obs)) for(name in ERGM_STATE_ELEMENTS) s.obs[[i]][[name]] <- state$s.obs.reduced[[i]][[name]]

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
      message("Saving state in ", sQuote(sprintf(control$checkpoint, iteration)),".")
      s.reduced <- s
      for(i in seq_along(s.reduced)) s.reduced[[i]]$model <- s.reduced[[i]]$proposal <- NULL
      if(obs){
        s.obs.reduced <- s.obs
        for(i in seq_along(s.obs.reduced)) s.obs.reduced[[i]]$model <- s.obs.reduced[[i]]$proposal <- NULL
      }
      save(list=intersect(ls(), c("s.reduced", "s.obs.reduced", STATE_VARIABLES, CONTROL_VARIABLES)), file=sprintf(control$checkpoint, iteration))
      rm(s.reduced)
      suppressWarnings(rm(s.obs.reduced))
    }

    # Obtain MCMC sample
    if(verbose) message("Starting unconstrained MCMC...")
    z <- ergm_MCMC_sample(s, control, theta=mcmc.init, verbose=max(verbose-1,0))

    if(z$status==1) stop("Number of edges in a simulated network exceeds that in the observed by a factor of more than ",floor(control$MCMLE.density.guard),". This is a strong indicator of model degeneracy or a very poor starting parameter configuration. If you are reasonably certain that neither of these is the case, increase the MCMLE.density.guard control.ergm() parameter.")
        
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero).
    statsmatrices <- z$stats
    s.returned <- z$networks
    statsmatrix <- as.matrix(statsmatrices)
    
    if(verbose){
      message("Back from unconstrained MCMC.")
      if(verbose>1){
        message("Average statistics:")
        message_print(colMeans(statsmatrix))
      }
    }

    ## Compute the sample estimating equations.
    esteqs <- ergm.estfun(statsmatrices, theta=mcmc.init, model=model)
    esteq <- as.matrix(esteqs)
    if(is.const.sample(esteq) && !all(esteq==0))
      stop("Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.")

    check_nonidentifiability(esteq, NULL, model,
                             tol = control$MCMLE.nonident.tol, type="statistics",
                             nonident_action = control$MCMLE.nonident,
                             nonvar_action = control$MCMLE.nonvar)

    ##  Do the same, if observation process:
    if(obs){
      if(adapt.obs.var) esteq.var <- var(esteq)
      .set_obs_samplesize()

      if(verbose) message("Starting constrained MCMC...")
      z.obs <- ergm_MCMC_sample(s.obs, control.obs, theta=mcmc.init, verbose=max(verbose-1,0))

      if(adapt.obs.var){
        obs.ESS.adj <- z.obs$final.effectiveSize / control$MCMC.effectiveSize
        if(verbose>1) message("New constrained vs. unconstrained target ESS adjustment factor: ", format(obs.ESS.adj), ".")
      }

      statsmatrices.obs <- z.obs$stats
      s.obs.returned <- z.obs$networks
      statsmatrix.obs <- as.matrix(statsmatrices.obs)
      
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
    
    if(control$MCMLE.sequential) {
      s <- s.returned
      
      if(obs){
        s.obs <- s.obs.returned
      }      
    }

    if(!is.null(control$MCMLE.save_intermediates)){
      save(list=intersect(ls(), INTERMEDIATE_VARIABLES), file=sprintf(control$MCMLE.save_intermediates, iteration))
    }

    esteqs.obs <- if(obs) ergm.estfun(statsmatrices.obs, theta=mcmc.init, model=model) else NULL
    esteq.obs <- if(obs) as.matrix(esteqs.obs) else NULL

    # Update the interval to be used.
    if(adapt){
      control$MCMC.interval <- round(max(z$final.interval/control$MCMLE.effectiveSize.interval_drop,1))
      control$MCMC.burnin <- round(max(z$final.interval*16,16))
      if(verbose) message("New interval = ",control$MCMC.interval,".")
      if(obs){
        control.obs$MCMC.interval <- round(max(z.obs$final.interval/control$MCMLE.effectiveSize.interval_drop,1))
        control.obs$MCMC.burnin <- round(max(z.obs$final.interval*16,16))
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
      s.returned <- lapply(s.returned, as.network)
      l <- list(coefficients=mcmc.init, mc.se=rep(NA,length=length(mcmc.init)),
                sample=statsmatrices, sample.obs=statsmatrices.obs,
                iterations=1, MCMCtheta=mcmc.init,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(mcmc.init)), #acf=NULL,
                samplesize=control$MCMC.samplesize, failure=TRUE,
                newnetwork = s.returned[[1]],
                newnetworks = s.returned)
      return(l)
    } 

    # Need to compute MCMC SE for "confidence" termination criterion
    # if it has the possibility of terminating.
    if(control$MCMLE.termination=='confidence'){
      estdiff <- NVL3(esteq.obs, colMeans(.), 0) - colMeans(esteq)
      pprec <- diag(sqrt(control$MCMLE.MCMC.precision), nrow=length(estdiff))
      Vm <- pprec%*%(cov(esteq) - NVL3(esteq.obs, cov(.), 0))%*%pprec
      novar <- diag(Vm) <= 0
      Vm[,novar] <- Vm[novar,] <- 0
      d2 <- xTAx_seigen(estdiff, Vm)
      if(d2<2) last.adequate <- TRUE
    }

      if(verbose){message("Starting MCMLE Optimization...")}

      if(!is.null(control$MCMLE.steplength.margin)){
        steplen <- .Hummel.steplength(
          if(control$MCMLE.steplength.esteq) esteq else statsmatrix[,!model$etamap$offsetmap,drop=FALSE],
          if(control$MCMLE.steplength.esteq) esteq.obs else statsmatrix.obs[,!model$etamap$offsetmap,drop=FALSE],
          control$MCMLE.steplength.margin, control$MCMLE.steplength, verbose=verbose,
          x2.num.max=control$MCMLE.steplength.miss.sample, parallel=control$MCMLE.steplength.parallel, control=control
        )

        # If the step length margin is negative and signals convergence,
        # rerun with margin of 0 and use the results to test
        # convergence.
        steplen0 <-
          if(control$MCMLE.termination%in%c("precision","Hummel") && control$MCMLE.steplength.margin<0 && control$MCMLE.steplength==steplen)
          .Hummel.steplength(
              if(control$MCMLE.steplength.esteq) esteq else statsmatrix[,!model$etamap$offsetmap,drop=FALSE],
              if(control$MCMLE.steplength.esteq) esteq.obs else statsmatrix.obs[,!model$etamap$offsetmap,drop=FALSE],
              0, control$MCMLE.steplength, verbose=verbose,
            x2.num.max=control$MCMLE.steplength.miss.sample,
            parallel=control$MCMLE.steplength.parallel, control=control
          )
          else steplen
        
        steplen.converged <- control$MCMLE.steplength==steplen0
      
      
      }else{
        steplen <- control$MCMLE.steplength
        steplen.converged <- TRUE
      }

      message("Optimizing with step length ", fixed.pval(steplen, eps = control$MCMLE.steplength.min), ".")
      if(control$MCMLE.steplength==steplen && !steplen.converged)
        message("Note that convergence diagnostic step length is ",steplen0,".")
      
        
      if(steplen.converged || is.null(control$MCMLE.steplength.margin) || iteration==control$MCMLE.maxit) calc.MCSE <- TRUE
      
      steplen.hist <- c(steplen.hist, steplen)
      
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(init=mcmc.init, model=model,
                       statsmatrices=statsmatrices, 
                       statsmatrices.obs=statsmatrices.obs, 
                       epsilon=control$epsilon,
                       nr.maxit=control$MCMLE.NR.maxit,
                       nr.reltol=control$MCMLE.NR.reltol,
                       calc.mcmc.se=control$MCMLE.termination == "precision" || (control$MCMC.addto.se && last.adequate) || iteration == control$MCMLE.maxit,
                       hessianflag=control$main.hessian,
                       method=control$MCMLE.method,
                       dampening=control$MCMLE.dampening,
                       dampening.min.ess=control$MCMLE.dampening.min.ess,
                       dampening.level=control$MCMLE.dampening.level,
                       metric=control$MCMLE.metric,
                       steplen=steplen,
                       verbose=verbose,
                       estimateonly=!calc.MCSE)
        message("The log-likelihood improved by ", fixed.pval(v$loglikelihood, 4), ".")
          
    coef.hist <- rbind(coef.hist, coef(v))
    stats.obs.hist <- NVL3(statsmatrix.obs, rbind(stats.obs.hist, apply(.[], 2, base::mean)))
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, base::mean))

    # This allows premature termination.
    
    if(control$MCMLE.termination=='Hotelling'){
      conv.pval <- ERRVL(try(suppressWarnings(approx.hotelling.diff.test(esteqs, esteqs.obs)$p.value)), NA)
      message("Nonconvergence test p-value:", format(conv.pval), "")
      # I.e., so that the probability of one false nonconvergence in two successive iterations is control$MCMLE.conv.min.pval (sort of).
      if(!is.na(conv.pval) && conv.pval>=1-sqrt(1-control$MCMLE.conv.min.pval)){
        if(last.adequate){
          message("No nonconvergence detected twice. Stopping.")
          break
        }else{
          message("No nonconvergence detected once; increasing sample size if not already increased.")
          last.adequate <- TRUE
          .boost_samplesize(control$MCMLE.last.boost, TRUE)
        }
      }else{
        last.adequate <- FALSE
      }
    }else if(control$MCMLE.termination=='confidence'){
      if(!is.null(estdiff.prev)){
        d2.prev <- xTAx_seigen(estdiff.prev, Vm)
        if(verbose) message("Distance from origin on tolerance region scale: ", format(d2), " (previously ", format(d2.prev), ").")
        d2.not.improved <- d2.not.improved[-1] 
        if(d2 >= d2.prev){
          d2.not.improved <- c(d2.not.improved,TRUE)
        }else{
          d2.not.improved <- c(d2.not.improved,FALSE)
        }
      }
      estdiff.prev <- estdiff
      
      if(d2<2){
        IS.lw <- function(sm, etadiff){
          nochg <- etadiff==0 | apply(sm, 2, function(x) max(x)==min(x)) # Also takes care of offset terms.
          # Calculate log-importance-weights
          basepred <- sm[,!nochg,drop=FALSE] %*% etadiff[!nochg]
        }
        lw2w <- function(lw){w<-exp(lw-max(lw)); w/sum(w)}
        # Handle a corner case in which none of the constrained sample
        # statistics vary. Then, Hotelling's test must be performed in
        # the 1-sample mode.
        novar.obs <- NVL3(esteq.obs, all(sweep(., 2L, .[1L,], `==`)), FALSE)
        hotel <- try(suppressWarnings(
          approx.hotelling.diff.test(esteqs,
                                     if(!novar.obs) esteqs.obs,
                                     mu0 = if(novar.obs) esteq.obs[1L,])
        ), silent=TRUE)
        if(inherits(hotel, "try-error")){ # Within tolerance ellipsoid, but cannot be tested.
          message("Unable to test for convergence; increasing sample size.")
          .boost_samplesize(control$MCMLE.confidence.boost)
        }else{
          etadiff <- ergm.eta(coef(v), model$etamap) - ergm.eta(mcmc.init, model$etamap)
          esteq.lw <- IS.lw(statsmatrix, etadiff)
          esteq.w <- lw2w(esteq.lw)
          estdiff <- -lweighted.mean(esteq, esteq.lw)
          estcov <- hotel$covariance.x*sum(esteq.w^2)*length(esteq.w)
          if(obs && !novar.obs){
            esteq.obs.lw <- IS.lw(statsmatrix.obs, etadiff)
            esteq.obs.w <- lw2w(esteq.obs.lw)
            estdiff <- estdiff + lweighted.mean(esteq.obs, esteq.obs.lw)
            estcov <- estcov + hotel$covariance.y*sum(esteq.obs.w^2)*length(esteq.obs.w)
          }
          estdiff <- estdiff[!hotel$novar]
          estcov <- estcov[!hotel$novar, !hotel$novar]

          d2e <- xTAx_seigen(estdiff, Vm[!hotel$novar, !hotel$novar])
          if(d2e<1){ # Update ends within tolerance ellipsoid.
            T2 <- try(.ellipsoid_mahalanobis(estdiff, estcov, Vm[!hotel$novar, !hotel$novar]), silent=TRUE) # Distance to the nearest point on the tolerance region boundary.
            if(inherits(T2, "try-error")){ # Within tolerance ellipsoid, but cannot be tested.
              message("Unable to test for convergence; increasing sample size.")
              .boost_samplesize(control$MCMLE.confidence.boost)
            }else{ # Within tolerance ellipsoid, can be tested.
              T2nullity <- attr(T2,"nullity")
              T2param <- hotel$parameter["param"] - T2nullity
              if(T2nullity && verbose) message("Estimated covariance matrix of the statistics has nullity ", format(T2nullity), ". Effective parameter number adjusted to ", format(T2param), ".")
              nonconv.pval <- .ptsq(T2, T2param, hotel$parameter["df"], lower.tail=FALSE)
              if(verbose) message("Test statistic: T^2 = ", format(T2),", with ",
                                  format(T2param), " free parameter(s) and ", format(hotel$parameter["df"]), " degrees of freedom.")
              message("Convergence test p-value: ", fixed.pval(nonconv.pval, 4), ". ", appendLF=FALSE)
              if(nonconv.pval < 1-control$MCMLE.confidence){
                message("Converged with ",control$MCMLE.confidence*100,"% confidence.")
                break
              }else{
                message("Not converged with ",control$MCMLE.confidence*100,"% confidence; increasing sample size.")
                critval <- .qtsq(control$MCMLE.confidence, T2param, hotel$parameter["df"])
                if(verbose) message(control$MCMLE.confidence*100,"% confidence critical value = ",format(critval),".")
                boost <- min((critval/T2),control$MCMLE.confidence.boost) # I.e., we want to increase the denominator far enough to reach the critical value.
                .boost_samplesize(boost)
              }
            }
          }
        }
      }
      # If either the estimating function is far from the tolerance
      # region *or* if it's close, but did not end the iteration
      # inside it.
      if(d2>=2 || d2e>1){ 
        message("Estimating equations are not within tolerance region.")
        if(sum(d2.not.improved) > control$MCMLE.confidence.boost.threshold){
          message("Estimating equations did not move closer to tolerance region more than ", control$MCMLE.confidence.boost.threshold," time(s) in ", control$MCMLE.confidence.boost.lag, " steps; increasing sample size.")
          .boost_samplesize(control$MCMLE.confidence.boost)
          d2.not.improved[] <- FALSE
        }
      }
    }else if(!steplen.converged){ # If step length is less than its maximum, don't bother with precision stuff.
      last.adequate <- FALSE
      .boost_samplesize(1, TRUE)
    }else if(control$MCMLE.termination == "precision"){
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
        
        if (adapt) { # ESS-based sampling
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
          if (adapt.obs.ESS) { # ESS-based sampling
            control.obs$MCMC.effectiveSize <- round(control.obs$MCMC.effectiveSize * prec.scl)
            if(control.obs$MCMC.effectiveSize/control.obs$MCMC.samplesize>control.obs$MCMLE.MCMC.max.ESS.frac) control.obs$MCMC.samplesize <- control.obs$MCMC.effectiveSize/control.obs$MCMLE.MCMC.max.ESS.frac
            # control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
            message("Increasing target constrained MCMC sample size to ", control.obs$MCMC.samplesize, ", ESS to",control.obs$MCMC.effectiveSize,".")
          } else if(adapt.obs.var){
            .set_obs_samplesize()
          } else { # Fixed-interval sampling
            control.obs$MCMC.samplesize <- round(control.obs$MCMC.samplesize * prec.scl)
            control.obs$MCMC.burnin <- round(control.obs$MCMC.burnin * prec.scl)
            message("Increasing constrained MCMC sample size to ", control.obs$MCMC.samplesize, ", burn-in to",control.obs$MCMC.burnin,".")
          }
        }
      }
    }else if(control$MCMLE.termination=='Hummel'){
      if(last.adequate){
        message("Step length converged twice. Stopping.")
        break
      }else{
        message("Step length converged once. Increasing MCMC sample size.")
        last.adequate <- TRUE
        .boost_samplesize(control$MCMLE.last.boost, TRUE)
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
    # Update the coefficient for MCMC sampling.
    mcmc.init <- coef(v)
  } # end of main loop

  message("Finished MCMLE.")

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- statsmatrices
  if(obs) v$sample.obs <- statsmatrices.obs
  nws.returned <- lapply(s.returned, as.network)
  v$newnetworks <- nws.returned
  v$newnetwork <- nws.returned[[1]]
  v$coef.init <- init
  v$est.cov <- v$mc.cov
  v$mc.cov <- NULL

  v$coef.hist <- coef.hist
  v$stats.hist <- stats.hist
  v$stats.obs.hist <- stats.obs.hist
  v$steplen.hist <- steplen.hist
  
  v$iterations <- iteration

  if(obs) for(name in OBS_MCMC_CONTROLS) control[[paste0("obs.", name)]] <- control.obs[[name]]
  v$control <- control
  
  v$etamap <- model$etamap
  v$MCMCflag <- TRUE
  v
}

#' Find the shortest squared Mahalanobis distance (with covariance W)
#' from a point `y` to an ellipsoid defined by `x'U^-1 x = 1`, provided
#' that `y` is in the interior of the ellipsoid.
#'
#' @param y a vector
#' @param W,U a square matrix
#'
#' @noRd
.ellipsoid_mahalanobis <- function(y, W, U, tol=sqrt(.Machine$double.eps)){
  y <- c(y)
  if(xTAx_seigen(y,U,tol=tol)>=1) stop("Point is not in the interior of the ellipsoid.")
  I <- diag(length(y))
  WUi <- t(qr.solve(qr(U, tol=tol), W))
  x <- function(l) c(solve(I+l*WUi, y)) # Singluar for negative reciprocals of eigenvalues of WiU.
  zerofn <- function(l) ERRVL(try({x <- x(l); xTAx_seigen(x,U,tol=tol)-1}, silent=TRUE), +Inf)

  # For some reason, WU sometimes has 0i element in its eigenvalues.
  eig <- Re(eigen(WUi, only.values=TRUE)$values)
  lmin <- -1/max(eig)
  l <- uniroot(zerofn, lower=lmin, upper=0, tol=sqrt(.Machine$double.xmin))$root
  x <- x(l)

  xTAx_seigen(y-x, W, tol=tol)
}
