############################################################################
# The <ergm.CD> function provides one of the styles of maximum
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
#                    this is ignored; default=control$CD.sequential
#   estimate       : whether to optimize the init coefficients via
#                    <ergm.estimate>; default=TRUE
#   ...            : additional parameters that may be passed from within;
#                    all are ignored
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.CD>=*);
#      note that if the model is degenerate, only 'coef' and 'sample' are
#      returned; if 'estimate'=FALSE, the MCMC and se variables will be
#      NA or NULL
#
#############################################################################

ergm.CD.fixed <- function(init, nw, model,
                             control, 
                             MHproposal, MHproposal.obs,
                             verbose=FALSE,
                             estimate=TRUE,
                             response=NULL, ...) {
  cat("Starting contrastive divergence estimation via CD-MCMLE:\n",sep="")
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  steplen.hist <- c()

  nthreads <- max(
    if(inherits(control$parallel,"cluster")) nrow(summary(control$parallel))
    else control$parallel,
    1)

  # Store information about original network, which will be returned at end
  nw.orig <- network.copy(nw)

  # Impute missing dyads.
  nw <- single.impute.dyads(nw, response=response)
  model$nw.stats <- summary(model$formula, response=response, basis=nw)

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
  if(obs){
    control$CD.nsteps<-control$CD.nsteps.obs
    control$CD.multiplicity<-control$CD.multiplicity.obs
  }
  
  # Initialize control.obs and other *.obs if there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$MCMC.samplesize <- control$obs.MCMC.samplesize
    control.obs$MCMC.interval <- control$obs.MCMC.interval
    control.obs$MCMC.burnin <- control$obs.MCMC.burnin
    control.obs$MCMC.burnin.min <- control$obs.MCMC.burnin.min

    nws.obs <- lapply(nws, network.copy)
    statshifts.obs <- statshifts
  }
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  finished <- FALSE

  for(iteration in 1:control$CD.maxit){
    if(iteration == control$CD.maxit) finished <- TRUE
    if(verbose){
      cat("Iteration ",iteration," of at most ", control$CD.maxit,
          " with parameter: \n", sep="")
      print(mcmc.init)
    }else{
      cat("Iteration ",iteration," of at most ", control$CD.maxit,": \n",sep="")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.init, model$etamap)
    z <- ergm.getCDsample(nws, model, MHproposal, mcmc.eta0, control, verbose, response=response, theta=mcmc.init, etamap=model$etamap)

    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrices <- mapply(sweep, z$statsmatrices, statshifts, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
    for(i in seq_along(statsmatrices)) colnames(statsmatrices[[i]]) <- model$coef.names
    statsmatrix <- do.call("rbind",statsmatrices)
    
    if(verbose){
      cat("Back from unconstrained CD. Average statistics:\n")
      print(apply(statsmatrix, 2, mean))
    }
    
    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm.getCDsample(nws.obs, model, MHproposal.obs, mcmc.eta0, control.obs, verbose, response=response, theta=mcmc.init, etamap=model$etamap)

      statsmatrices.obs <- mapply(sweep, z.obs$statsmatrices, statshifts.obs, MoreArgs=list(MARGIN=2, FUN="+"), SIMPLIFY=FALSE)
      for(i in seq_along(statsmatrices.obs)) colnames(statsmatrices.obs[[i]]) <- model$coef.names
      statsmatrix.obs <- do.call("rbind",statsmatrices.obs)
      
      if(verbose){
        cat("Back from constrained MCMC. Average statistics:\n")
        print(apply(statsmatrix.obs, 2, mean))
      }
    }else{
      statsmatrices.obs <- statsmatrix.obs <- NULL
      z.obs <- NULL
    }

    # Compute the sample estimating equations and the convergence p-value. 
    esteq <- .ergm.esteq(mcmc.init, model, statsmatrix)
    if(isTRUE(all.equal(apply(esteq,2,sd), rep(0,ncol(esteq)), check.names=FALSE))&&!all(esteq==0))
      stop("Unconstrained CD sampling did not mix at all. Optimization cannot continue.")
    esteq.obs <- if(obs) .ergm.esteq(mcmc.init, model, statsmatrix.obs) else NULL   
    conv.pval <- suppressWarnings(approx.hotelling.diff.test(esteq, esteq.obs, assume.indep=TRUE)$p.value)
                                            
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..
    if(verbose){
      cat("Average estimating equation values:\n")
      print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else colMeans(esteq))
    }
    cat("Convergence test P-value:",format(conv.pval, scientific=TRUE,digits=2),"\n")
    if(conv.pval>control$CD.conv.min.pval){
      cat("Convergence detected. Stopping.\n")
      finished <- TRUE
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
                newnetwork = nw)
      return(structure (l, class="ergm"))
    } 

    statsmatrix.0 <- statsmatrix
    statsmatrix.0.obs <- statsmatrix.obs
    if(control$CD.steplength=="adaptive"){
      if(verbose){cat("Calling adaptive MCMLE Optimization...\n")}
      adaptive.steplength <- 2
      statsmean <- apply(statsmatrix.0,2,mean)
      v <- list(loglikelihood=control$CD.adaptive.trustregion*2)
      while(v$loglikelihood > control$CD.adaptive.trustregion){
        adaptive.steplength <- adaptive.steplength / 2
        if(!is.null(statsmatrix.0.obs)){
          statsmatrix.obs <- sweep(statsmatrix.0.obs,2,(colMeans(statsmatrix.0.obs)-statsmean)*(1-adaptive.steplength))
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
                         nr.maxit=control$CD.NR.maxit,
                         nr.reltol=control$CD.NR.reltol,
                         calc.mcmc.se=FALSE, hessianflag=control$main.hessian,
                         trustregion=control$CD.trustregion, method=control$CD.method,
                         metric=control$CD.metric,
                         dampening=control$CD.dampening,
                         dampening.min.ess=control$CD.dampening.min.ess,
                         dampening.level=control$CD.dampening.level,
                         compress=control$MCMC.compress, verbose=verbose,
                         estimateonly=TRUE)
      }
      if(v$loglikelihood < control$CD.trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        cat("The log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("The log-likelihood did not improve.\n")
      }
      steplen.hist <- c(steplen.hist, adaptive.steplength)
    }else{
      steplen <-
        if(!is.null(control$CD.steplength.margin))
          .Hummel.steplength(
            if(control$MCMLE.Hummel.esteq) esteq else statsmatrix.0[,!model$etamap$offsetmap,drop=FALSE], 
            if(control$MCMLE.Hummel.esteq) esteq.obs else statsmatrix.0.obs[,!model$etamap$offsetmap,drop=FALSE],
            control$CD.steplength.margin, control$CD.steplength)
        else control$CD.steplength
      
      if(verbose){cat("Calling MCMLE Optimization...\n")}
      statsmean <- apply(statsmatrix.0,2,mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- sweep(statsmatrix.0.obs,2,(colMeans(statsmatrix.0.obs)-statsmean)*(1-steplen))
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
                       nr.maxit=control$CD.NR.maxit,
                       nr.reltol=control$CD.NR.reltol,
                       calc.mcmc.se=FALSE, 
                       hessianflag=control$main.hessian,
                       trustregion=control$CD.trustregion, 
                       method=control$CD.method,
                       dampening=control$CD.dampening,
                       dampening.min.ess=control$CD.dampening.min.ess,
                       dampening.level=control$CD.dampening.level,
                       metric=control$CD.metric,
                       compress=control$MCMC.compress, verbose=verbose,
                       estimateonly=!finished)
      if(v$loglikelihood < control$CD.trustregion-0.001){
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
    stats.obs.hist <- if(!is.null(statsmatrix.obs)) rbind(stats.obs.hist, apply(statsmatrix.obs[], 2, mean)) else NULL
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, mean))
    if(finished) break # This allows premature termination.
  } # end of main loop

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, control)
  
  v$network <- nw.orig
  v$newnetwork <- nw
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

