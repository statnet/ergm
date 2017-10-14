#  File R/ergm.MCMLE.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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
#   model          : the model, as returned by <ergm.getmodel>
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
#   MHproposal     : an MHproposal object for 'nw', as returned by
#                    <MHproposal>
#   MHproposal.obs : an MHproposal object for the observed network of'nw',
#                    as returned by <MHproposal>
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
                             MHproposal, MHproposal.obs,
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
  control$obs.MCMC.effectiveSize <- control$obs.MCMLE.effectiveSize
  
  control$MCMC.base.effectiveSize <- control$MCMC.effectiveSize
  control$obs.MCMC.base.effectiveSize <- control$obs.MCMC.effectiveSize
  
  control$MCMC.base.samplesize <- control$MCMC.samplesize
  control$obs.MCMC.base.samplesize <- control$obs.MCMC.samplesize

  nthreads <- max(
    if(inherits(control$parallel,"cluster")) nrow(summary(control$parallel))
    else control$parallel,
    1)
  
  # Store information about original network, which will be returned at end
  nw.orig <- network.copy(nw)

  # Impute missing dyads.
  nw <- single.impute.dyads(nw, response=response)
  model$nw.stats <- ergm.getglobalstats(nw, model, response=response)

  if(control$MCMLE.density.guard>1){
    # Calculate the density guard threshold.
    control$MCMC.max.maxedges <- round(min(control$MCMC.max.maxedges,
                                           max(control$MCMLE.density.guard*network.edgecount(nw,FALSE),
                                               control$MCMLE.density.guard.min)))
    control$MCMC.init.maxedges <- round(min(control$MCMC.max.maxedges, control$MCMC.init.maxedges))
    if(verbose) message("Density guard set to ",control$MCMC.max.maxedges," from an initial count of ",network.edgecount(nw,FALSE)," edges.")
  }  

  nws <- rep(list(nw),nthreads) # nws is now a list of networks.

  # statshift is the difference between the target.stats (if
  # specified) and the statistics of the networks in the LHS of the
  # formula or produced by SAN. If target.stats is not speficied
  # explicitly, they are computed from this network, so
  # statshift==0. To make target.stats play nicely with offsets, we
  # set statshifts to 0 where target.stats is NA (due to offset).
  statshift <- model$nw.stats - NVL(model$target.stats,model$nw.stats)
  statshift[is.na(statshift)] <- 0
  statshifts <- rep(list(statshift), nthreads) # Each network needs its own statshift.

  # Is there observational structure?
  obs <- ! is.null(MHproposal.obs)
  
  # Initialize control.obs and other *.obs if there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$MCMC.base.samplesize <- control$obs.MCMC.base.samplesize
    control.obs$MCMC.base.effectiveSize <- control$obs.MCMC.base.effectiveSize
    control.obs$MCMC.samplesize <- control$obs.MCMC.samplesize
    control.obs$MCMC.effectiveSize <- control$obs.MCMC.effectiveSize
    control.obs$MCMC.interval <- control$obs.MCMC.interval
    control.obs$MCMC.burnin <- control$obs.MCMC.burnin

    nws.obs <- lapply(nws, network::network.copy)
    statshifts.obs <- statshifts
  }

  # A helper function to increase the MCMC sample size and target effective size by the specified factor.
  .boost_samplesize <- function(boost, base=FALSE){
    control <- get("control", parent.frame())
    control$MCMC.samplesize <- round((if(base) control$MCMC.base.samplesize else control$MCMC.samplesize) * boost)
    control$MCMC.effectiveSize <- NVL3((if(base) control$MCMC.base.effectiveSize else control$MCMC.effectiveSize), . * boost)
    assign("control", control, parent.frame())
    if(obs){
      control.obs <- get("control.obs", parent.frame())
      control.obs$MCMC.samplesize <- round((if(base) control.obs$MCMC.base.samplesize else control.obs$MCMC.samplesize) * boost)
      control.obs$MCMC.effectiveSize <- NVL3((if(base) control.obs$MCMC.base.effectiveSize else control.obs$MCMC.effectiveSize), . * boost)
      assign("control.obs", control.obs, parent.frame())
    }
    NULL
  }


  if(control$MCMLE.termination=='confidence'){
    estdiff.prev <- NULL
  }

  
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  calc.MCSE <- FALSE
  last.adequate <- FALSE
  
  for(iteration in 1:control$MCMLE.maxit){
    if(verbose){
      message("\nIteration ",iteration," of at most ", control$MCMLE.maxit,
          " with parameter:")
      message_print(mcmc.init)
    }else{
      message("Iteration ",iteration," of at most ", control$MCMLE.maxit,":")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.init, model$etamap)
    z <- ergm.getMCMCsample(nws, model, MHproposal, mcmc.eta0, control, verbose, response=response, theta=mcmc.init, etamap=model$etamap, update.nws=FALSE)
        
    if(z$status==1) stop("Number of edges in a simulated network exceeds that in the observed by a factor of more than ",floor(control$MCMLE.density.guard),". This is a strong indicator of model degeneracy or a very poor starting parameter configuration. If you are reasonably certain that neither of these is the case, increase the MCMLE.density.guard control.ergm() parameter.")
        
    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrices <- mapply(sweep, z$statsmatrices, statshifts, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
    for(i in seq_along(statsmatrices)) colnames(statsmatrices[[i]]) <- model$coef.names
    nws.returned <- z$newnetworks
    statsmatrix <- do.call(rbind,statsmatrices)
    
    if(verbose){
      message("Back from unconstrained MCMC. Average statistics:")
      message_print(apply(statsmatrix, 2, base::mean))
    }
    
    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm.getMCMCsample(nws.obs, NVL(model$obs.model,model), MHproposal.obs, mcmc.eta0, control.obs, verbose, response=response, theta=mcmc.init, etamap=model$etamap, update.nws=FALSE)
      
      if(z.obs$status==1) stop("Number of edges in the simulated network exceeds that observed by a large factor (",control$MCMC.max.maxedges,"). This is a strong indication of model degeneracy. If you are reasonably certain that this is not the case, increase the MCMLE.density.guard control.ergm() parameter.")
      
      statsmatrices.obs <- mapply(sweep, z.obs$statsmatrices, statshifts.obs, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
      for(i in seq_along(statsmatrices.obs)) colnames(statsmatrices.obs[[i]]) <- model$coef.names
      nws.obs.returned <- z.obs$newnetworks
      statsmatrix.obs <- do.call(rbind,statsmatrices.obs)
      
      if(verbose){
        message("Back from constrained MCMC. Average statistics:")
        message_print(apply(statsmatrix.obs, 2, base::mean))
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

    # Compute the sample estimating equations and the convergence p-value. 
    esteq <- .ergm.esteq(mcmc.init, model, statsmatrix)
    if(isTRUE(all.equal(apply(esteq,2,stats::sd), rep(0,ncol(esteq)), check.names=FALSE))&&!all(esteq==0))
      stop("Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.")
    esteq.obs <- if(obs) .ergm.esteq(mcmc.init, model, statsmatrix.obs) else NULL

    # Update the interval to be used.
    if(!is.null(control$MCMC.effectiveSize)){
      control$MCMC.interval <- round(max(z$final.interval,2)/2)
      if(verbose) message("New interval = ",control$MCMC.interval,".")
      if(obs){
        control.obs$MCMC.interval <- round(max(z.obs$final.interval,2)/2)
        if(verbose) message("New constrained interval = ",control.obs$MCMC.interval,".")
      }
    }
        
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..

    # These are only nontrivial when the model is curved or when there are missing data.
    if(verbose && (is.curved(model)||obs)){
      message("Average estimating equation values:")
      message_print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else colMeans(esteq))
    }

    if(!estimate){
      if(verbose){message("Skipping optimization routines...")}
      nws.returned <- lapply(nws.returned, newnw.extract, response=response)
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


    # Need to compute MCMC SE for "confidence" termination criterion
    # if it has the possibility of terminating.
    if(control$MCMLE.termination=='confidence'){
      estdiff <- NVL3(esteq.obs, colMeans(.), 0) - colMeans(esteq)
      pprec <- diag(sqrt(control$MCMLE.MCMC.precision), nrow=length(estdiff))
      Vm <- pprec%*%(cov(esteq) - NVL3(esteq.obs, cov(.), 0))%*%pprec
      Vm <- as.matrix(nearPD(Vm, posd.tol=0)$mat) # Ensure tolerance hyperellipsoid is PSD. (If it's not PD, the ellipsoid is workable, if flat.)
      iVm <- ginv(Vm)
      d2 <- estdiff%*%iVm%*%estdiff
      if(d2<1) last.adequate <- TRUE
    }

    if(control$MCMLE.steplength=="adaptive"){
      if(verbose){message("Starting adaptive MCMLE Optimization...")}
      adaptive.steplength <- 2
      v <- list(loglikelihood=control$MCMLE.adaptive.trustregion*2)
      while(v$loglikelihood > control$MCMLE.adaptive.trustregion){
        adaptive.steplength <- adaptive.steplength / 2
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
                         steplen=adaptive.steplength,
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
      
      if(!is.null(control$MCMLE.steplength.margin)){
        steplen <- .Hummel.steplength(
          if(control$MCMLE.Hummel.esteq) esteq else statsmatrix[,!model$etamap$offsetmap,drop=FALSE], 
          if(control$MCMLE.Hummel.esteq) esteq.obs else statsmatrix.obs[,!model$etamap$offsetmap,drop=FALSE],
          control$MCMLE.steplength.margin, control$MCMLE.steplength,point.gamma.exp=control$MCMLE.steplength.point.exp,steplength.prev=steplen,x1.prefilter=control$MCMLE.steplength.prefilter,x2.prefilter=control$MCMLE.steplength.prefilter,verbose=verbose,
          x2.num.max=control$MCMLE.Hummel.miss.sample, steplength.maxit=control$MCMLE.Hummel.maxit,
          last=(iteration==control$MCMLE.maxit))

        # If the step length margin is negative and signals convergence,
        # rerun with margin of 0 and use the results to test
        # convergence.
        steplen0 <-
          if(control$MCMLE.termination%in%c("precision","Hummel") && control$MCMLE.steplength.margin<0 && control$MCMLE.steplength==steplen)
            .Hummel.steplength(
              if(control$MCMLE.Hummel.esteq) esteq else statsmatrix[,!model$etamap$offsetmap,drop=FALSE], 
              if(control$MCMLE.Hummel.esteq) esteq.obs else statsmatrix.obs[,!model$etamap$offsetmap,drop=FALSE],
              0, control$MCMLE.steplength,steplength.prev=steplen,point.gamma.exp=control$MCMLE.steplength.point.exp,x1.prefilter=control$MCMLE.steplength.prefilter,x2.prefilter=control$MCMLE.steplength.prefilter,verbose=verbose,
              x2.num.max=control$MCMLE.Hummel.miss.sample, steplength.maxit=control$MCMLE.Hummel.maxit,
              last=(iteration==control$MCMLE.maxit))
          else steplen
        
        steplen.converged <- control$MCMLE.steplength==steplen0
        
        
      }else{
        steplen <- control$MCMLE.steplength
        steplen.converged <- TRUE
      }

      message("Optimizing with step length ",steplen,".")
      if(control$MCMLE.steplength==steplen && !steplen.converged)
        message("Note that convergence diagnostic step length is ",steplen0,".")
      
        
      if(steplen.converged || is.null(control$MCMLE.steplength.margin) || iteration==control$MCMLE.maxit) calc.MCSE <- TRUE
      
      steplen.hist <- c(steplen.hist, steplen)
      
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
                       steplen=steplen, steplen.point.exp=control$MCMLE.steplength.point.exp,
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
    
    if(control$MCMLE.termination=='Hotelling'){
      conv.pval <- ERRVL(try(approx.hotelling.diff.test(esteq, esteq.obs)$p.value), NA)
      message("Nonconvergence test p-value:",conv.pval,"")
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
      if(is.null(iVm)){
        message("Tolerance region could not be computed; increasing sample size.")
        .boost_samplesize(control$MCMLE.confidence.boost)
      }else{
        if(d2>=1){ # Not within tolerance ellipsoid.
          message("Estimating equations are not within tolerance region.")
          if(!is.null(estdiff.prev)){
            d2.prev <- estdiff.prev%*%iVm%*%estdiff.prev
            if(verbose) message("Distance from origin on tolerance region scale: ", d2, " (previously ", d2.prev, ").")
            if(d2 > d2.prev){
              message("Estimating equations did not move closer to tolerance region; increasing sample size.")
              .boost_samplesize(control$MCMLE.confidence.boost)
            }
          }
        }else{
          hotel <- try(approx.hotelling.diff.test(esteq, esteq.obs))
          if(inherits(hotel, "try-error")){ # Within tolerance ellipsoid, but cannot be tested.
            message("Unable to test for convergence; increasing sample size.")
            .boost_samplesize(control$MCMLE.confidence.boost)
          }else{ # Within tolerance ellipsoid, can be tested.
            T2 <- with(hotel, .ellipsoid_mahalanobis(estimate, covariance, iVm[!novar, !novar])) # Distance to the nearest point on the tolerance region boundary.
            nonconv.pval <- .ptsq(T2, hotel$parameter["param"], hotel$parameter["df"], lower.tail=FALSE)
            if(verbose) message("Test statistic: T^2 = ",T2,", with ",
                                hotel$parameter["param"], " free parameters and ",hotel$parameter["df"], " degrees of freedom.")
            message("Convergence test p-value: ",nonconv.pval,". ", appendLF=FALSE)
            if(nonconv.pval < 1-control$MCMLE.confidence){
              message("Converged with ",control$MCMLE.confidence*100,"% confidence.")
              break
            }else{
              message("Not converged with ",control$MCMLE.confidence*100,"% confidence; increasing sample size.")
              critval <- .qtsq(control$MCMLE.confidence, hotel$parameter["param"], hotel$parameter["df"])
              if(verbose) message(control$MCMLE.confidence*100,"% confidence critical value = ",critval,".")
              boost <- min((critval/T2),control$MCMLE.confidence.boost) # I.e., we want to increase the denominator far enough to reach the critical value.
              .boost_samplesize(boost)
            }
          }
        }
        estdiff.prev <- estdiff
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
          if (!is.null(control.obs$MCMC.effectiveSize)) { # ESS-based sampling
            control.obs$MCMC.effectiveSize <- round(control.obs$MCMC.effectiveSize * prec.scl)
            if(control.obs$MCMC.effectiveSize/control.obs$MCMC.samplesize>control.obs$MCMLE.MCMC.max.ESS.frac) control.obs$MCMC.samplesize <- control.obs$MCMC.effectiveSize/control.obs$MCMLE.MCMC.max.ESS.frac
            # control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
            message("Increasing target constrained MCMC sample size to ", control.obs$MCMC.samplesize, ", ESS to",control.obs$MCMC.effectiveSize,".")
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
  } # end of main loop

  message("Finished MCMLE.")

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.obs, control)
  
  nws.returned <- lapply(nws.returned, newnw.extract, response=response)
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

#' Find the shortest squared Mahalanobis distance (with covariance W)
#' from a point `y` to an ellipsoid defined by `x'U x = 1`, provided
#' that `y` is in the interior of the ellipsoid.
#'
#' @param y a vector
#' @param W,U a square matrix
#'
#' @noRd
.ellipsoid_mahalanobis <- function(y, W, U){
  y <- c(y)
  if(y%*%U%*%y>=1) stop("Point is not in the interior of the ellipsoid.")
  I <- diag(length(y))
  WU <- W%*%U
  x <- function(l) c(solve(I+l*WU, y)) # Singluar for negative reciprocals of eigenvalues of WiU.
  zerofn <- function(l) {x <- x(l); c(x%*%U%*%x)-1}

  # For some reason, WU sometimes has 0i element in its eigenvalues.
  eig <- Re(eigen(WU, only.values=TRUE)$values)
  lmin <- -1/max(eig)+sqrt(.Machine$double.eps)
  l <- uniroot(zerofn, lower=lmin, upper=0)$root
  x <- x(l)
  (y-x)%*%solve(W)%*%(y-x)
}
