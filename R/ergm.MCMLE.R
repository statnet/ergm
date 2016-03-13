#  File R/ergm.MCMLE.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
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
  cat("Starting maximum likelihood estimation via MCMLE:\n",sep="")
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  steplen.hist <- c()
  steplen <- control$MCMLE.steplength
  if(control$MCMLE.steplength=="adaptive") steplen <- 1

  control$MCMC.effectiveSize <- control$MCMLE.effectiveSize
  control$MCMC.base.samplesize <- control$MCMC.samplesize

  nthreads <- max(
    if(inherits(control$parallel,"cluster")) nrow(summary(control$parallel))
    else control$parallel,
    1)
  
  # Store information about original network, which will be returned at end
  nw.orig <- network.copy(nw)

  # Impute missing dyads.
  nw <- single.impute.dyads(nw, response=response)
  model$nw.stats <- summary(model$formula, response=response, basis=nw)

  if(control$MCMLE.density.guard>1){
    # Calculate the density guard threshold.
    control$MCMC.max.maxedges <- round(min(control$MCMC.max.maxedges,
                                           max(control$MCMLE.density.guard*network.edgecount(nw,FALSE),
                                               control$MCMLE.density.guard.min)))
    control$MCMC.init.maxedges <- round(min(control$MCMC.max.maxedges, control$MCMC.init.maxedges))
    if(verbose) cat("Density guard set to",control$MCMC.max.maxedges,"from an initial count of",network.edgecount(nw,FALSE)," edges.\n")
  }  

  nws <- rep(list(nw),nthreads) # nws is now a list of networks.

  # statshift is the difference between the target.stats (if
  # specified) and the statistics of the networks in the LHS of the
  # formula or produced by SAN. If target.stats is not speficied
  # explicitly, they are computed from this network, so
  # statshift==0. To make target.stats play nicely with offsets, we
  # set statshifts to 0 where target.stats is NA (due to offset).
  statshift <- model$nw.stats - model$target.stats
  statshift[is.na(statshift)] <- 0
  statshifts <- rep(list(statshift), nthreads) # Each network needs its own statshift.

  # Is there observational structure?
  obs <- ! is.null(MHproposal.obs)
  
  # Initialize control.obs and other *.obs if there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$MCMC.samplesize <- control$obs.MCMC.samplesize
    control.obs$MCMC.interval <- control$obs.MCMC.interval
    control.obs$MCMC.burnin <- control$obs.MCMC.burnin
    control.obs$MCMC.burnin.min <- control$obs.MCMC.burnin.min

    nws.obs <- lapply(nws, network::network.copy)
    statshifts.obs <- statshifts
  }
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  calc.MCSE <- FALSE
  last.adequate <- FALSE
  
  for(iteration in 1:control$MCMLE.maxit){
    if(verbose){
      cat("Iteration ",iteration," of at most ", control$MCMLE.maxit,
          " with parameter: \n", sep="")
      print(mcmc.init)
    }else{
      cat("Iteration ",iteration," of at most ", control$MCMLE.maxit,": \n",sep="")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.init, model$etamap)
    z <- ergm.getMCMCsample(nws, model, MHproposal, mcmc.eta0, control, verbose, response=response, theta=mcmc.init, etamap=model$etamap)
        
    if(z$status==1) stop("Number of edges in a simulated network exceeds that in the observed by a factor of more than ",floor(control$MCMLE.density.guard),". This is a strong indicator of model degeneracy or a very poor starting parameter configuration. If you are reasonably certain that neither of these is the case, increase the MCMLE.density.guard control.ergm() parameter.")
        
    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrices <- mapply(sweep, z$statsmatrices, statshifts, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
    for(i in seq_along(statsmatrices)) colnames(statsmatrices[[i]]) <- model$coef.names
    nws.returned <- lapply(z$newnetworks,network::network.copy)
    statsmatrix <- do.call(rbind,statsmatrices)
    
    if(verbose){
      cat("Back from unconstrained MCMC. Average statistics:\n")
      print(apply(statsmatrix, 2, base::mean))
    }
    
    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm.getMCMCsample(nws.obs, model, MHproposal.obs, mcmc.eta0, control.obs, verbose, response=response, theta=mcmc.init, etamap=model$etamap)
      
      if(z.obs$status==1) stop("Number of edges in the simulated network exceeds that observed by a large factor (",control$MCMC.max.maxedges,"). This is a strong indication of model degeneracy. If you are reasonably certain that this is not the case, increase the MCMLE.density.guard control.ergm() parameter.")
      
      statsmatrices.obs <- mapply(sweep, z.obs$statsmatrices, statshifts.obs, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
      for(i in seq_along(statsmatrices.obs)) colnames(statsmatrices.obs[[i]]) <- model$coef.names
      nws.obs.returned <- lapply(z.obs$newnetworks, network::network.copy)
      statsmatrix.obs <- do.call(rbind,statsmatrices.obs)
      
      if(verbose){
        cat("Back from constrained MCMC. Average statistics:\n")
        print(apply(statsmatrix.obs, 2, base::mean))
      }
    }else{
      statsmatrices.obs <- statsmatrix.obs <- NULL
      z.obs <- NULL
    }
    
    if(sequential) {
      nws <- nws.returned
      statshifts <- lapply(nws, function(nw) summary(model$formula, basis=nw, response=response) - model$target.stats)
      
      if(obs){
        nws.obs <- nws.obs.returned
        statshifts.obs <- lapply(nws.obs, function(nw.obs) summary(model$formula, basis=nw.obs, response=response) - model$target.stats)
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
      if(verbose) cat("New interval =",control$MCMC.interval,".\n")
      if(obs){
        control.obs$MCMC.interval <- round(max(z.obs$final.interval,2)/2)
        if(verbose) cat("New constrained interval =",control.obs$MCMC.interval,".\n")
      }
    }
        
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..
    if(verbose){
      cat("Average estimating equation values:\n")
      print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else colMeans(esteq))
    }

    if(!estimate){
      if(verbose){cat("Skipping optimization routines...\n")}
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
      if(verbose){cat("Calling adaptive MCMLE Optimization...\n")}
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
        if(verbose){cat(paste("Using Newton-Raphson Step with step length",adaptive.steplength,"...\n"))}
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
        cat("The log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("The log-likelihood did not improve.\n")
      }
      steplen.hist <- c(steplen.hist, adaptive.steplength)
      steplen <- adaptive.steplength
    }else{
      steplen <-
        if(!is.null(control$MCMLE.steplength.margin))
          .Hummel.steplength(
            if(control$MCMLE.Hummel.esteq) esteq else statsmatrix.0[,!model$etamap$offsetmap,drop=FALSE], 
            if(control$MCMLE.Hummel.esteq) esteq.obs else statsmatrix.0.obs[,!model$etamap$offsetmap,drop=FALSE],
            control$MCMLE.steplength.margin, control$MCMLE.steplength,steplength.prev=steplen,verbose=verbose,
            x2.num.max=control$MCMLE.Hummel.miss.sample, steplen.maxit=control$MCMLE.Hummel.maxit)
        else control$MCMLE.steplength
      
      if(steplen==control$MCMLE.steplength || is.null(control$MCMLE.steplength.margin) || iteration==control$MCMLE.maxit) calc.MCSE <- TRUE
      
      if(verbose){cat("Calling MCMLE Optimization...\n")}
      statsmean <- apply(statsmatrix.0,2,base::mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- t(steplen*t(statsmatrix.0.obs) + (1-steplen)*statsmean) # I.e., shrink each point of statsmatrix.obs towards the centroid of statsmatrix.
      }else{
        statsmatrix <- sweep(statsmatrix.0,2,(1-steplen)*statsmean,"-")
      }
      steplen.hist <- c(steplen.hist, steplen)
      
      if(verbose){cat(paste("Using Newton-Raphson Step with step length ",steplen," ...\n"))}
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
        cat("The log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("The log-likelihood did not improve.\n")
      }
    }
          
    mcmc.init <- v$coef
    coef.hist <- rbind(coef.hist, mcmc.init)
    stats.obs.hist <- if(!is.null(statsmatrix.obs)) rbind(stats.obs.hist, apply(statsmatrix.obs[], 2, base::mean)) else NULL
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, base::mean))
    
    # This allows premature termination.
    
    if(steplen<control$MCMLE.steplength){ # If step length is less than its maximum, don't bother with precision stuff.
      last.adequate <- FALSE
      control$MCMC.samplesize <- control$MCMC.base.samplesize
      
    } else {
    
    if(control$MCMLE.termination == "precision"){
      prec.loss <- (sqrt(diag(v$mc.cov+v$covar))-sqrt(diag(v$covar)))/sqrt(diag(v$mc.cov+v$covar))
      if(verbose){
        cat("Standard Error:\n")
        print(sqrt(diag(v$covar)))
        cat("MC SE:\n")
        print(sqrt(diag(v$mc.cov)))
        cat("Linear scale precision loss due to MC estimation of the likelihood:\n")
        print(prec.loss)
      }
      if(sqrt(mean(prec.loss^2, na.rm=TRUE)) <= control$MCMLE.MCMC.precision){
        if(last.adequate){
          cat("Precision adequate twice. Stopping.\n")
          break
        }else{
          cat("Precision adequate. Performing one more iteration.\n")
          last.adequate <- TRUE
        }
      }else{
        last.adequate <- FALSE
        prec.scl <- max(sqrt(mean(prec.loss^2, na.rm=TRUE))/control$MCMLE.MCMC.precision, 1) # Never decrease it.
        if (!is.null(control$MCMC.effectiveSize)) { # ESS-based sampling
          control$MCMC.effectiveSize <- round(control$MCMC.effectiveSize * prec.scl)
          if(control$MCMC.effectiveSize/control$MCMC.samplesize>control$MCMLE.MCMC.max.ESS.frac) control$MCMC.samplesize <- control$MCMC.effectiveSize/control$MCMLE.MCMC.max.ESS.frac
          # control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
          cat("Increasing target MCMC sample size to ", control$MCMC.samplesize, ", ESS to",control$MCMC.effectiveSize,".\n")
        } else { # Fixed-interval sampling
          control$MCMC.samplesize <- round(control$MCMC.samplesize * prec.scl)
          control$MCMC.burnin <- round(control$MCMC.burnin * prec.scl)
          cat("Increasing MCMC sample size to ", control$MCMC.samplesize, ", burn-in to",control$MCMC.burnin,".\n")
        }
      }
    }else if(control$MCMLE.termination=='Hotelling'){
      conv.pval <- ERRVL(try(approx.hotelling.diff.test(esteq, esteq.obs)$p.value), NA)
      cat("Nonconvergence test p-value:",conv.pval,"\n")
      if(!is.na(conv.pval) && conv.pval>=control$MCMLE.conv.min.pval){
        cat("No nonconvergence detected. Stopping.\n")
        break
      }      
    }else if(control$MCMLE.termination=='Hummel'){
      if(last.adequate){
        cat("Step length converged twice. Stopping.\n")
        break
      }else{
        cat("Step length converged once. Increasing MCMC sample size.\n")
        last.adequate <- TRUE
        control$MCMC.samplesize <- control$MCMC.base.samplesize * control$MCMLE.last.boost
      }
    }
    
    }
    
    # stop if MCMLE is stuck (steplen stuck near 0)
    if ((length(steplen.hist) > 2) && sum(tail(steplen.hist,2)) < 2*control$MCMLE.steplength.min) {
      stop("MCMLE estimation stuck. There may be excessive correlation between model terms, suggesting a poor model for the observed data. If target.stats are specified, try increasing SAN parameters.")
    }    
    #Otherwise, don't stop before iterations are exhausted.
    if (iteration == control$MCMLE.maxit) {
      message("MCMLE estimation did not converge after ", control$MCMLE.maxit, " iterations. The estimated coefficients may not be accurate. Estimation may be resumed by passing the coefficients as initial values; see 'init' under ?control.ergm for details.\n")
    }
  } # end of main loop

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, control)
  
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

