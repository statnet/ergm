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

ergm.CD <- function(init, nw, model,
                             control, 
                             MHproposal, MHproposal.obs,
                             verbose=FALSE,
                             sequential=control$MCMLE.sequential,
                             estimate=TRUE,
                             response=NULL, ...) {
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  tether.hist <- 0
  
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

  # statshift is the difference between the target.stats (if
  # specified) and the statistics of the networks in the LHS of the
  # formula or produced by SAN. If target.stats is not speficied
  # explicitly, they are computed from this network, so statshift==0.
  statshift <- model$nw.stats - model$target.stats

  # Is there observational structure?
  obs <- ! is.null(MHproposal.obs)
  
  # Initialize control.obs in case there is observation structure
  
  if(obs){
    control.obs <- control
    control.obs$MCMC.samplesize <- control$obs.MCMC.samplesize
    control.obs$MCMC.interval <- control$obs.MCMC.interval
    control.obs$MCMC.burnin <- control$obs.MCMC.burnin

    nw.obs <- network.copy(nw)
    statshift.obs <- statshift
  }
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  parametervalues <- init # Keep track of all parameter values

  iteration <- 1
  repeat{ # Tether length loop: finish when doubling tether length doesn't change statistics
    cat("Tether length ", control$CD.nsteps,".\n", sep="")
    MCMLE.starting <- TRUE
    MCMLE.converged <- FALSE
  repeat{ # Optimization loop: finish when simulated is statistically indistinguishable from observed
    finished <- FALSE
    if(iteration == control$CD.maxit) finished <- TRUE

    tether.hist <- c(tether.hist, control$CD.nsteps)
    
    if(verbose){
      cat("Optimization run ",iteration," of at most ", control$CD.maxit, " with tether length ", control$CD.nsteps, " and with parameter: \n", sep="")
      print(mcmc.init)
    }else{
      cat("Optimization run ",iteration," of at most ", control$CD.maxit,": \n",sep="")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.init, model$etamap)
    if(control$CD.nsteps==Inf){
    z <- ergm.getMCMCsample(nw, model, MHproposal, mcmc.eta0, control, verbose, response=response, theta=mcmc.init, etamap=model$etamap)
    if(z$status==1) stop("Number of edges in a simulated network exceeds that in the observed by a factor of more than ",floor(control$MCMLE.density.guard),". This is a strong indicator of model degeneracy. If you are reasonably certain that this is not the case, increase the MCMLE.density.guard control.ergm() parameter.")
    }
    else
      z <- ergm.getCDsample(nw, model, MHproposal, mcmc.eta0, control, verbose, response=response)
    
    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
    colnames(statsmatrix) <- model$coef.names
    
    if(control$CD.nsteps==Inf)
    nw.returned <- network.copy(z$newnetwork)

    if(verbose){
      cat("Back from unconstrained MCMC. Average statistics:\n")
      print(apply(statsmatrix, 2, mean))
    }
    
    ##  Does the same, if observation process:
    if(obs){
      if(control$CD.nsteps==Inf){
      z.obs <- ergm.getMCMCsample(nw.obs, model, MHproposal.obs, mcmc.eta0, control.obs, verbose, response=response, theta=mcmc.init, etamap=model$etamap)

      if(z.obs$status==1) stop("Number of edges in the simulated network exceeds that observed by a large factor (",control$MCMC.max.maxedges,"). This is a strong indication of model degeneracy. If you are reasonably certain that this is not the case, increase the MCMLE.density.guard control.ergm() parameter.")
      }else
        z.obs <- ergm.getCDsample(nw.obs, model, MHproposal.obs, mcmc.eta0, control.obs, verbose, response=response)
      
      statsmatrix.obs <- sweep(z.obs$statsmatrix, 2, statshift.obs, "+")
      colnames(statsmatrix.obs) <- model$coef.names
      if(control$CD.nsteps==Inf)
      nw.obs.returned <- network.copy(z.obs$newnetwork)
      
      if(verbose){
        cat("Back from constrained MCMC. Average statistics:\n")
        print(apply(statsmatrix.obs, 2, mean))
      }
    }else{
      statsmatrix.obs <- NULL
    }
    
    if(control$CD.nsteps==Inf && sequential) {
      nw <- nw.returned
      statshift <- summary(model$formula, basis=nw, response=response) - model$target.stats

      if(obs){
        nw.obs <- nw.obs.returned
        statshift.obs <- summary(model$formula, basis=nw.obs, response=response) - model$target.stats
      }      
    }

    if(MCMLE.starting && control$CD.nsteps>1){
      if(obs){
        # There is certainly a better way to do this:
        smu <- unique(statsmatrix)
        smou <- unique(statsmatrix.obs)
        if(!any(apply(smu, 1, is.inCH, smou))){
          cat("Convex hulls of constrained and unconstrained sample statistics do not overlap. Reducing tether length.")
          control$CD.nsteps <- round(control$CD.nsteps * 3/4)
          break
        }       
      }else{
        if(!is.inCH(rep(0,ncol(statsmatrix)),statsmatrix)){
          cat("Convex hull of the sample does not contain the observed statistics. Reducing tether length.")
          control$CD.nsteps <- round(control$CD.nsteps * 3/4)
          break
        }
      }
    }

    # Compute the sample estimating equations and the convergence p-value.
    esteq <- .ergm.esteq(mcmc.init, model, statsmatrix)
    esteq.obs <- if(obs) .ergm.esteq(mcmc.init, model, statsmatrix.obs) else NULL   
    conv.pval <- approx.hotelling.diff.test(esteq, esteq.obs)$p.value
    
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..
    
    if(verbose){
      cat("Average estimating equation values:\n")
        print(if(obs) colMeans(esteq.obs)-colMeans(esteq) else colMeans(esteq))
    }
    cat("Convergence test P-value:",format(conv.pval, scientific=TRUE,digits=2),"\n")
    if(conv.pval>control$CD.conv.min.pval){
      MCMLE.converged <- TRUE
      if(control$CD.nsteps==Inf){
        cat("Convergence detected for untethered MCMC. Stopping.\n")
        finished <- TRUE
      }else if(MCMLE.starting && control$CD.nsteps>=control$CD.min.nsteps){
        cat("Convergence in tether lengths detected. Proceeding to untethered MCMC.\n")
        control$CD.nsteps <- Inf
      }else{
        cat("Convergence for this tether length detected. Increasing tether length.\n")
        control$CD.nsteps <- control$CD.nsteps*2
        if(control$CD.nsteps>=control$MCMC.interval/2){
          cat("Tether length exceeds the half the MCMC interval setting. Switching to untethered MCMC.\n")
          control$CD.nsteps <- Inf
        }
      }
    }
    MCMLE.starting <- FALSE

    if(!estimate){
      if(verbose){cat("Skipping optimization routines...\n")}
      l <- list(coef=mcmc.init, mc.se=rep(NA,length=length(mcmc.init)),
                sample=statsmatrix, sample.obs=statsmatrix.obs,
                iterations=1, MCMCtheta=mcmc.init,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(mcmc.init)), #acf=NULL,
                samplesize=control$MCMC.samplesize, failure=TRUE,
                newnetwork = nw.returned)
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
          statsmatrix.obs <- statsmatrix.0.obs*adaptive.steplength+statsmatrix.0*(1-adaptive.steplength)
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
                         calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
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
    }else{
      
      if(verbose){cat("Calling MCMLE Optimization...\n")}
      statsmean <- apply(statsmatrix.0,2,mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- statsmatrix.0.obs*control$CD.steplength+statsmatrix.0*(1-control$CD.steplength)
      }else{
        statsmatrix <- sweep(statsmatrix.0,2,(1-control$CD.steplength)*statsmean,"-")
      }
      if(verbose){cat(paste("Using Newton-Raphson Step with step length ",control$CD.steplength," ...\n"))}
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(init=mcmc.init, model=model,
                       statsmatrix=statsmatrix, 
                       statsmatrix.obs=statsmatrix.obs, 
                       epsilon=control$epsilon,
                       nr.maxit=control$CD.NR.maxit,
                       nr.reltol=control$CD.NR.reltol,
                       calc.mcmc.se=control$MCMC.addto.se, 
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
    stats.obs.hist <- if(!is.null(statsmatrix.obs)) rbind(stats.obs.hist, apply(statsmatrix.obs, 2, mean)) else NULL
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, mean))
    parametervalues <- rbind(parametervalues, mcmc.init)

    iteration <- iteration + 1
    
    if(MCMLE.converged || finished) break
  } # End of the optimization loop
    
    if(finished) break
  } # End of the step length loop

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, control)
  
  v$network <- nw.orig
  v$newnetwork <- if(control$CD.nsteps==Inf) nw.returned else nw.orig
  v$coef.init <- init
  
  v$coef.hist <- coef.hist
  v$stats.hist <- stats.hist
  v$stats.obs.hist <- stats.obs.hist
  v$tether.hist <- tether.hist
  # The following output is sometimes helpful.  It's the total history
  # of all eta values, from the initial eta0 to the final estimate
  # v$allparamvals <- parametervalues


  v$etamap <- model$etamap
  v
}

