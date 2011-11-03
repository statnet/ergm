############################################################################
# The <ergm.mainfitloop> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is the default and uses
# optimization of an MCMC estimate of the log-likelihood.  (The other
# MLE styles are found in functions <ergm.robmon>, <ergm.stocapprox>, and
# <ergm.stepping> 
#
#
# --PARAMETERS--
#   theta0         : the initial theta values
#   nw             : the network 
#   model          : the model, as returned by <ergm.getmodel>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
#   initialfit     : an ergm object, as the initial fit, possibly returned
#                    by <ergm.initialfit>
#   MCMCparams     : a list of parameters for controlling the MCMC sampling;
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
#                    this is ignored; default=MCMCparams$sequential
#   estimate       : whether to optimize the theta0 coefficients via
#                    <ergm.estimate>; default=TRUE
#   ...            : additional parameters that may be passed from within;
#                    all are ignored
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.mainfitloop>=*);
#      note that if the model is degenerate, only 'coef' and 'sample' are
#      returned; if 'estimate'=FALSE, the MCMC and se variables will be
#      NA or NULL
#
#############################################################################

ergm.mainfitloop <- function(theta0, nw, model, Clist,
                             initialfit, 
                             MCMCparams, 
                             MHproposal, MHproposal.obs,
                             verbose=FALSE,
                             sequential=MCMCparams$sequential,
                             estimate=TRUE,
                             response=NULL, ...) {
  # Initialize the history of parameters and statistics.
  theta.hist <- rbind(theta0)
  stats.hist <- matrix(NA, 0, length(Clist$obs))
  stats.obs.hist <- matrix(NA, 0, length(Clist$obs))
  
  # Store information about original network, which will be returned at end
  nw.orig <- network.copy(nw)

  # Calculate the amount by which all of the MCMC statistics should be adjusted
  # to account for the fact that they are all calculated relative to the
  # observed network.  Unless otherwise specified, the value of meanstats
  # is simply the observed statistics, which means statshift equals zero
  # most of the time.
  statshift <- Clist$obs - Clist$meanstats
  MCMCparams$meanstats <- Clist$meanstats

  # Is there observational structure?
  obs <- ! is.null(MHproposal.obs)
  
  # Initialize MCMCparams.obs in case there is observation structure
  
  if(obs){
    MCMCparams.obs <- MCMCparams
    if(!is.null(MCMCparams$obs.MCMCsamplesize)){
      MCMCparams.obs$MCMCsamplesize <- MCMCparams$obs.MCMCsamplesize
    }
    if(!is.null(MCMCparams$obs.interval)){
      MCMCparams.obs$interval <- MCMCparams$obs.interval
    }
    if(!is.null(MCMCparams$obs.burnin)){
      MCMCparams.obs$burnin <- MCMCparams$obs.burnin
    }
  }
  iteration <- 0
  finished <- FALSE
  # mcmc.theta0 will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  theta0 will never change.
  mcmc.theta0 <- theta0
  parametervalues <- theta0 # Keep track of all parameter values
  while(!finished){
	  iteration <- iteration + 1
    if(verbose){
      cat("Iteration ",iteration," of at most ", MCMCparams$maxit,
          " with parameter: \n", sep="")
      print(mcmc.theta0)
    }else{
      cat("Iteration ",iteration," of at most ", MCMCparams$maxit,": \n",sep="")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.theta0, model$etamap)
    z <- ergm.getMCMCsample.parallel(nw, model, MHproposal, mcmc.eta0, MCMCparams, verbose, response=response)
    
    # post-processing of sample statistics:  Shift each row by the
    # vector Clist$obs - Clist$meanstats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative meanstats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
    colnames(statsmatrix) <- model$coef.names
    nw.returned <- network.copy(z$newnetwork)

    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm.getMCMCsample.parallel(nw, model, MHproposal.obs, mcmc.eta0, MCMCparams.obs, verbose, response=response)
      statsmatrix.obs <- sweep(z.obs$statsmatrix, 2, statshift, "+")
      colnames(statsmatrix.obs) <- model$coef.names
      nw.obs.returned <- network.copy(z.obs$newnetwork)
      if(verbose){cat("Back from constrained MCMC...\n")}
    }else{
      statsmatrix.obs <- NULL
      if(verbose){cat("Back from unconstrained MCMC...\n")}
      if(sequential) {
        nw <- nw.returned
        nw.obs <- summary(model$formula, basis=nw, response=response)
        namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
        statshift <- -Clist$meanstats
        statshift[!is.na(namesmatch)] <- statshift[!is.na(namesmatch)] + nw.obs[namesmatch[!is.na(namesmatch)]]
      }
    }
    
    # Removed block A of code here.  (See end of file.)
    
    if(!estimate){
      if(verbose){cat("Skipping optimization routines...\n")}
      l <- list(coef=mcmc.theta0, mc.se=rep(NA,length=length(mcmc.theta0)),
                sample=statsmatrix, sample.obs=statsmatrix.obs,
                iterations=1, MCMCtheta=mcmc.theta0,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(mcmc.theta0)), #acf=NULL,
                samplesize=MCMCparams$samplesize, failure=TRUE,
                newnetwork = nw.returned)
      return(structure (l, class="ergm"))
    } 

    statsmatrix.0 <- statsmatrix
    statsmatrix.0.obs <- statsmatrix.obs
    if(MCMCparams$steplength=="adaptive"){
      if(verbose){cat("Calling adaptive MCMLE Optimization...\n")}
      adaptive.steplength <- 2
      statsmean <- apply(statsmatrix.0,2,mean)
      v <- list(loglikelihood=MCMCparams$adaptive.trustregion*2)
      while(v$loglikelihood > MCMCparams$adaptive.trustregion){
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
        v<-ergm.estimate(theta0=mcmc.theta0, model=model,
                         statsmatrix=statsmatrix, 
                         statsmatrix.obs=statsmatrix.obs, 
                         epsilon=MCMCparams$epsilon,
                         nr.maxit=MCMCparams$nr.maxit,
                         nr.reltol=MCMCparams$nr.reltol,
                         calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                         trustregion=MCMCparams$trustregion, method=MCMCparams$method,
                         metric=MCMCparams$metric,
                         compress=MCMCparams$compress, verbose=verbose,
                         estimateonly=TRUE)
      }
      if(v$loglikelihood < MCMCparams$trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        cat("the log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("the log-likelihood did not improve.\n")
      }
      if((adaptive.steplength==1) && (v$loglikelihood < MCMCparams$adaptive.epsilon) ){break}
    }else{

      if(verbose){cat("Calling MCMLE Optimization...\n")}
      statsmean <- apply(statsmatrix.0,2,mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- statsmatrix.0.obs*MCMCparams$steplength+statsmatrix.0*(1-MCMCparams$steplength)
      }else{
        statsmatrix <- sweep(statsmatrix.0,2,(1-MCMCparams$steplength)*statsmean,"-")
      }
      if(verbose){cat(paste("Using Newton-Raphson Step with step length ",MCMCparams$steplength," ...\n"))}
      finished <- iteration >= MCMCparams$maxit
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(theta0=mcmc.theta0, model=model,
                       statsmatrix=statsmatrix, 
                       statsmatrix.obs=statsmatrix.obs, 
                       epsilon=MCMCparams$epsilon,
                       nr.maxit=MCMCparams$nr.maxit,
                       nr.reltol=MCMCparams$nr.reltol,
                       calc.mcmc.se=MCMCparams$calc.mcmc.se, 
                       hessianflag=MCMCparams$hessian,
                       trustregion=MCMCparams$trustregion, 
                       method=MCMCparams$method,
                       metric=MCMCparams$metric,
                       compress=MCMCparams$compress, verbose=verbose,
                       estimateonly=!finished)
      if(v$loglikelihood < MCMCparams$trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        cat("the log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("the log-likelihood did not improve.\n")
      }
      if((MCMCparams$steplength==1) && (v$loglikelihood < MCMCparams$adaptive.epsilon) ){break}
    }
          
    finished <- iteration >= MCMCparams$maxit
    mcmc.theta0 <- v$coef
    theta.hist <- rbind(theta.hist, mcmc.theta0)
    stats.obs.hist <- if(!is.null(statsmatrix.obs)) rbind(stats.obs.hist, apply(statsmatrix.obs, 2, mean)) else NULL
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, mean))
    parametervalues <- rbind(parametervalues, mcmc.theta0)
  } # end of main loop

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, MCMCparams) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, MCMCparams)
  
  v$burnin <- MCMCparams$burnin
  v$samplesize <- MCMCparams$samplesize
  v$interval <- MCMCparams$interval
  v$network <- nw.orig
  v$newnetwork <- nw.returned
  v$interval <- MCMCparams$interval
  v$theta.original <- theta0
  v$mplefit <- initialfit
  v$parallel <- MCMCparams$parallel
  
  v$theta.hist <- theta.hist
  v$stats.hist <- stats.hist
  v$stats.obs.hist <- stats.obs.hist
  # The following output is sometimes helpful.  It's the total history
  # of all eta values, from the initial eta0 to the final estimate
  # v$allparamvals <- parametervalues


  v$null.deviance <- 2*network.dyadcount(nw.orig)*log(2)
  v$etamap <- model$etamap
  v
}

#################
# Block of code A removed from above:

##  Check for degeneracy if new network has fewer than 49999 edges
#    if(z$nedges >= 50000-1 || ergm.checkdegeneracy(statsmatrix, statsmatrix.obs, verbose=verbose)){
#      if(iteration <= MCMCparams$maxit){
#        cat(paste("The MCMC sampler is producing degenerate samples.\n",
#                  "Try starting the algorithm at an alternative model\n",
#                  "(That is, changing the 'theta0' argument).\n",
#                  "I am trying something simple...\n",
#                  "The current theta0 is:\n"))
#        print(mcmc.theta0)
#        mcmc.theta0 <- 0.9*mcmc.theta0
#        next
#      }else{
#        cat(paste("The MCMC sampler is producing degenerate samples.\n",
#                  "Try starting the algorithm at an alternative model\n",
#                  "(That is, changing the 'theta0' argument).\n",
#                  "The current theta0 is:\n"))
#        print(mcmc.theta0)
#        v$coef <- mcmc.theta0
#        return(structure (v, class="ergm"))
#      }
#    }
#    if(verbose){
#      cat(paste("The density of the returned network is",
#                network.density(nw.returned),"\n"))
#      cat(paste("The density of the original network is",
#                network.density(nw.orig),"\n"))
#      cat("Summary of simulation, relative to observed network:\n")
#      print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
#      degreedist(nw.returned)
#      cat("Meanstats of simulation, relative to observed network:\n")
#      print(summary(model$formula, basis=nw.returned)-Clist$meanstats)
#      if(network.naedgecount(nw) > 0){
#        cat("Summary of simulation, relative to missing network:\n")
#        a = apply(statsmatrix.miss,2,summary.statsmatrix.ergm)[4,]
#        b = sweep(apply(statsmatrix,2,summary.statsmatrix.ergm),2,a,"-")
#        print(b,scipen=6)
#        degreedist(nw.miss.returned)
#        cat("Meanstats of simulation, relative to missing network:\n")
#        print(summary(model$formula, basis=nw.miss.returned)-Clist$meanstats)
#        nw.returned <- network.copy(nw.miss.returned)
#      }
#    }


