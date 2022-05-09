#  File R/ergm.CD.fixed.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
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
                             proposal, proposal.obs,
                             verbose=FALSE,
                             estimate=TRUE, ...) {
  message("Starting contrastive divergence estimation via CD-MCMLE:")
  # Is there observational structure?
  obs <- ! is.null(proposal.obs)
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  steplen.hist <- c()
  steplen <- control$CD.steplength

  if(is.null(control$CD.samplesize)) control$CD.samplesize <- control$CD.samplesize.per_theta*nparam(model,canonical=FALSE, offset=FALSE)
  if(obs && is.null(control$obs.CD.samplesize)) control$obs.CD.samplesize <- control$obs.CD.samplesize.per_theta*nparam(model,canonical=FALSE, offset=FALSE)

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
  s <- single.impute.dyads(nw, constraints=proposal$arguments$constraints, constraints.obs=proposal.obs$arguments$constraints, min_informative = control$obs.MCMC.impute.min_informative, default_density = control$obs.MCMC.impute.default_density, output="ergm_state", verbose=verbose)

  # statshift is the difference between the target.stats (if
  # specified) and the statistics of the networks in the LHS of the
  # formula or produced by SAN. If target.stats is not speficied
  # explicitly, they are computed from this network, so
  # statshift==0. To make target.stats play nicely with offsets, we
  # set statshifts to 0 where target.stats is NA (due to offset).
  model$nw.stats <- summary(model, s)
  statshift <- model$nw.stats - NVL(model$target.stats,model$nw.stats)
  statshift[is.na(statshift)] <- 0
  s <- update(s, model=model, proposal=proposal, stats=statshift)

  s <- rep(list(s),nthreads(control)) # s is now a list of states.
  
  # Initialize control.obs and other *.obs if there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$CD.nsteps<-control$CD.nsteps.obs
    control.obs$CD.multiplicity<-control$CD.multiplicity.obs
    control.obs$CD.samplesize <- control$obs.CD.samplesize
    control.obs$CD.interval <- control$obs.CD.interval
    control.obs$CD.burnin <- control$obs.CD.burnin

    s.obs <- lapply(s, update, model=NVL(model$obs.model,model), proposal=proposal.obs)
  }
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the CD samples.  init will never change.
  mcmc.init <- init
  finished <- FALSE

  for(iteration in 1:control$CD.maxit){
    if(iteration == control$CD.maxit) finished <- TRUE
    if(verbose){
      message("\nIteration ",iteration," of at most ", control$CD.maxit,
          " with parameter:")
      message_print(mcmc.init)
    }else{
      message("Iteration ",iteration," of at most ", control$CD.maxit,":")
    }

    # Obtain CD sample
    z <- ergm_CD_sample(s, control, verbose=verbose, theta=mcmc.init)

    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero).
    statsmatrices <- z$stats
    statsmatrix <- as.matrix(statsmatrices)
    
    if(verbose){
      message("Back from unconstrained CD.")
      if(verbose>1){
        message("Average statistics:")
        message_print(colMeans(statsmatrix))
      }
    }
    
    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm_CD_sample(s.obs, control.obs, theta=mcmc.init, verbose=max(verbose-1,0))

      statsmatrices.obs <- z.obs$stats
      statsmatrix.obs <- as.matrix(statsmatrices.obs)
      
      if(verbose){
        message("Back from constrained CD.")
        if(verbose>1){
          message("Average statistics:")
          message_print(colMeans(statsmatrix.obs))
        }
      }
    }else{
      statsmatrices.obs <- statsmatrix.obs <- NULL
      z.obs <- NULL
    }

    # Compute the sample estimating functions and the convergence p-value. 
    esteqs <- ergm.estfun(statsmatrices, theta=mcmc.init, model=model)
    esteq <- as.matrix(esteqs)
    if(isTRUE(all.equal(apply(esteq,2,stats::sd), rep(0,ncol(esteq)), check.names=FALSE))&&!all(esteq==0))
      stop("Unconstrained CD sampling did not mix at all. Optimization cannot continue.")
    esteqs.obs <- if(obs) ergm.estfun(statsmatrices.obs, theta=mcmc.init, model=model) else NULL
    esteq.obs <- if(obs) as.matrix(esteqs.obs) else NULL

    conv.pval <- suppressWarnings(approx.hotelling.diff.test(esteq, esteq.obs, assume.indep=TRUE)$p.value)

    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..
    if(verbose){
      message("Average estimating function values:")
      message_print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else -colMeans(esteq))
    }
    message("Convergence test P-value:",format(conv.pval, scientific=TRUE,digits=2),"")
    if(conv.pval>control$CD.conv.min.pval){
      message("Convergence detected. Stopping.")
      finished <- TRUE
    }

    if(!estimate){
      if(verbose){message("Skipping optimization routines...")}
      l <- list(coefficients=mcmc.init, mc.se=rep(NA,length=length(mcmc.init)),
                sample=statsmatrices, sample.obs=statsmatrices.obs,
                iterations=1, MCMCtheta=mcmc.init,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(mcmc.init)), #acf=NULL,
                samplesize=control$CD.samplesize, failure=TRUE,
                newnetwork = nw,
                newnetworks = nw)
      return(structure (l, class="ergm"))
    } 

      if(verbose){message("Calling CD-MCMLE Optimization...")}
      steplen <-
        if(!is.null(control$CD.steplength.margin))
          .Hummel.steplength(
            if(control$CD.steplength.esteq) esteq else statsmatrix[,!model$etamap$offsetmap,drop=FALSE], 
            if(control$CD.steplength.esteq) esteq.obs else statsmatrix.obs[,!model$etamap$offsetmap,drop=FALSE],
            control$CD.steplength.margin, control$CD.steplength, verbose=verbose,
            x2.num.max=control$CD.steplength.miss.sample,
            parallel=control$CD.steplength.parallel, control=control)
        else control$CD.steplength
      
      steplen.hist <- c(steplen.hist, steplen)
      # stop if MCMLE is stuck (steplen stuck near 0)
      if ((length(steplen.hist) > 2) && sum(tail(steplen.hist,2)) < 2*control$CD.steplength.min) {
        stop("CD-MCMLE estimation stuck. There may be excessive correlation between model terms, suggesting a poor model for the observed data. If target.stats are specified, try increasing SAN parameters.")
      }    
      
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(init=mcmc.init, model=model,
                       statsmatrices=statsmatrices, 
                       statsmatrices.obs=statsmatrices.obs, 
                       epsilon=control$epsilon,
                       nr.maxit=control$CD.NR.maxit,
                       nr.reltol=control$CD.NR.reltol,
                       calc.mcmc.se=FALSE,
                       hessianflag=control$main.hessian,
                       method=control$CD.method,
                       dampening=control$CD.dampening,
                       dampening.min.ess=control$CD.dampening.min.ess,
                       dampening.level=control$CD.dampening.level,
                       metric=control$CD.metric,
                       steplen=steplen,
                       verbose=verbose,
                       estimateonly=!finished)
        current.scipen <- options()$scipen
        options(scipen=3)
        message("The log-likelihood improved by ",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),".")
        options(scipen=current.scipen)
          
    coef.hist <- rbind(coef.hist, coef(v))
    stats.obs.hist <- NVL3(statsmatrix.obs, rbind(stats.obs.hist, apply(.[], 2, base::mean)))
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, base::mean))
    if(finished) break # This allows premature termination.
    # Update the coefficient for CD sampling.
    mcmc.init <- coef(v)
  } # end of main loop

  message("Finished CD.")

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- statsmatrices
  if(obs) v$sample.obs <- statsmatrices.obs
  
  v$network <- nw.orig
  v$newnetworks <- nw
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

