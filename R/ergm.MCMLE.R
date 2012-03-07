#  File ergm/R/ergm.MCMLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
############################################################################
# The <ergm.MCMLE> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is the default and uses
# optimization of an MCMC estimate of the log-likelihood.  (The other
# MLE styles are found in functions <ergm.robmon>, <ergm.stocapprox>, and
# <ergm.stepping> 
#############################################################################

ergm.MCMLE <- function(init, nw, model,
                             initialfit, 
                             control, 
                             MHproposal, MHproposal.obs,
                             verbose=FALSE,
                             sequential=control$MCMLE.sequential,
                             estimate=TRUE, ...) {
  # Initialize the history of parameters and statistics.
  coef.hist <- rbind(init)
  stats.hist <- matrix(NA, 0, length(model$nw.stats))
  stats.obs.hist <- matrix(NA, 0, length(model$nw.stats))
  
  # Store information about original network, which will be returned at end
  nw.orig <- network.copy(nw)

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
    control.obs$MCMC.samplesize <- control$MCMLE.obs.MCMC.samplesize
    control.obs$MCMC.interval <- control$MCMLE.obs.MCMC.interval
    control.obs$MCMC.burnin <- control$MCMLE.obs.MCMC.burnin
  }
  finished <- FALSE
  # mcmc.init will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  init will never change.
  mcmc.init <- init
  parametervalues <- init # Keep track of all parameter values
  for(iteration in 1:control$MCMLE.maxit){
    if(iteration == control$MCMLE.maxit) finished <- TRUE
    if(verbose){
      cat("Iteration ",iteration," of at most ", control$MCMLE.maxit,
          " with parameter: \n", sep="")
      print(mcmc.init)
    }else{
      cat("Iteration ",iteration," of at most ", control$MCMLE.maxit,": \n",sep="")
    }

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.init, model$etamap)
    z <- ergm.getMCMCsample(nw, model, MHproposal, mcmc.eta0, control, verbose)
    
    # post-processing of sample statistics:  Shift each row by the
    # vector model$nw.stats - model$target.stats, store returned nw
    # The statistics in statsmatrix should all be relative to either the
    # observed statistics or, if given, the alternative target.stats
    # (i.e., the estimation goal is to use the statsmatrix to find 
    # parameters that will give a mean vector of zero)
    statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
    colnames(statsmatrix) <- model$coef.names
    nw.returned <- network.copy(z$newnetwork)

    if(verbose){
      cat("Back from unconstrained MCMC. Average statistics:\n")
      print(apply(statsmatrix, 2, mean))
    }
    
    ##  Does the same, if observation process:
    if(obs){
      z.obs <- ergm.getMCMCsample(nw, model, MHproposal.obs, mcmc.eta0, control.obs, verbose)
      statsmatrix.obs <- sweep(z.obs$statsmatrix, 2, statshift, "+")
      colnames(statsmatrix.obs) <- model$coef.names
      nw.obs.returned <- network.copy(z.obs$newnetwork)
      
      if(verbose){
        cat("Back from constrained MCMC. Average statistics:\n")
        print(apply(statsmatrix.obs, 2, mean))
      }
    }else{
      statsmatrix.obs <- NULL
      if(sequential) {
        nw <- nw.returned
        nw.obs <- summary(model$formula, basis=nw)
        namesmatch <- match(names(model$target.stats), names(nw.obs))
        statshift <- -model$target.stats
        statshift[!is.na(namesmatch)] <- statshift[!is.na(namesmatch)] + nw.obs[namesmatch[!is.na(namesmatch)]]
      }
    }

    # If the model is linear, all non-offset statistics are passed. If
    # the model is curved, the (likelihood) estimating equations (3.1)
    # by Hunter and Handcock (2006) are given instead.
    conv.pval <- approx.hotelling.diff.test(t(ergm.etagradmult(mcmc.init,t(statsmatrix),model$etamap))[,!model$etamap$offsettheta,drop=FALSE],

                                            if(obs) t(ergm.etagradmult(mcmc.init,t(statsmatrix.obs),model$etamap))[,!model$etamap$offsettheta,drop=FALSE])
    # We can either pretty-print the p-value here, or we can print the
    # full thing. What the latter gives us is a nice "progress report"
    # on whether the estimation is getting better..

    # if(verbose) cat("P-value for equality of observed and simulated statistics:",format.pval(conv.pval,digits=4,eps=1e-4),"\n")
    if(verbose) cat("P-value for equality of observed and simulated statistics:",conv.pval,"\n")
    if(conv.pval>control$MCMLE.conv.min.pval){
      cat("Convergence detected. Stopping early.\n")
      finished <- TRUE
    }
    
    # Removed block A of code here.  (See end of file.)
    
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
    if(control$MCMLE.steplength=="adaptive"){
      if(verbose){cat("Calling adaptive MCMLE Optimization...\n")}
      adaptive.steplength <- 2
      statsmean <- apply(statsmatrix.0,2,mean)
      v <- list(loglikelihood=control$MCMLE.adaptive.trustregion*2)
      while(v$loglikelihood > control$MCMLE.adaptive.trustregion){
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
                         nr.maxit=control$MCMLE.NR.maxit,
                         nr.reltol=control$MCMLE.NR.reltol,
                         calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
                         trustregion=control$MCMLE.trustregion, method=control$MCMLE.method,
                         metric=control$MCMLE.metric,
                         compress=control$MCMC.compress, verbose=verbose,
                         estimateonly=TRUE)
      }
      if(v$loglikelihood < control$MCMLE.trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        cat("the log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("the log-likelihood did not improve.\n")
      }
      if((adaptive.steplength==1) && (v$loglikelihood < control$MCMLE.adaptive.epsilon) ){break}
    }else{

      if(verbose){cat("Calling MCMLE Optimization...\n")}
      statsmean <- apply(statsmatrix.0,2,mean)
      if(!is.null(statsmatrix.0.obs)){
        statsmatrix.obs <- statsmatrix.0.obs*control$MCMLE.steplength+statsmatrix.0*(1-control$MCMLE.steplength)
      }else{
        statsmatrix <- sweep(statsmatrix.0,2,(1-control$MCMLE.steplength)*statsmean,"-")
      }
      if(verbose){cat(paste("Using Newton-Raphson Step with step length ",control$MCMLE.steplength," ...\n"))}
      # Use estimateonly=TRUE if this is not the last iteration.
      v<-ergm.estimate(init=mcmc.init, model=model,
                       statsmatrix=statsmatrix, 
                       statsmatrix.obs=statsmatrix.obs, 
                       epsilon=control$epsilon,
                       nr.maxit=control$MCMLE.NR.maxit,
                       nr.reltol=control$MCMLE.NR.reltol,
                       calc.mcmc.se=control$MCMC.addto.se, 
                       hessianflag=control$main.hessian,
                       trustregion=control$MCMLE.trustregion, 
                       method=control$MCMLE.method,
                       metric=control$MCMLE.metric,
                       compress=control$MCMC.compress, verbose=verbose,
                       estimateonly=!finished)
      if(v$loglikelihood < control$MCMLE.trustregion-0.001){
        current.scipen <- options()$scipen
        options(scipen=3)
        cat("the log-likelihood improved by",
            format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
        options(scipen=current.scipen)
      }else{
        cat("the log-likelihood did not improve.\n")
      }
      if((control$MCMLE.steplength==1) && (v$loglikelihood < control$MCMLE.adaptive.epsilon) ){break}
    }
          
    mcmc.init <- v$coef
    coef.hist <- rbind(coef.hist, mcmc.init)
    stats.obs.hist <- if(!is.null(statsmatrix.obs)) rbind(stats.obs.hist, apply(statsmatrix.obs, 2, mean)) else NULL
    stats.hist <- rbind(stats.hist, apply(statsmatrix, 2, mean))
    parametervalues <- rbind(parametervalues, mcmc.init)
    if(finished) break # This allows premature termination.
  } # end of main loop

  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- ergm.sample.tomcmc(statsmatrix.0, control) 
  if(obs) v$sample.obs <- ergm.sample.tomcmc(statsmatrix.0.obs, control)
  
  v$network <- nw.orig
  v$newnetwork <- nw.returned
  v$coef.init <- init
  v$initialfit <- initialfit
  
  v$coef.hist <- coef.hist
  v$stats.hist <- stats.hist
  v$stats.obs.hist <- stats.obs.hist
  # The following output is sometimes helpful.  It's the total history
  # of all eta values, from the initial eta0 to the final estimate
  # v$allparamvals <- parametervalues


  v$null.deviance <- 2*network.dyadcount(nw.orig)*log(2)
  v$etamap <- model$etamap
  v
}

approx.hotelling.diff.test<-function(x,y=NULL){
  d <-  colMeans(x)
  if(!is.null(y)) d <- d - colMeans(y)

  # If a statistic doesn't vary and doesn't match, don't terminate.
  x.n <- effectiveSize(x)
  if(any(d[x.n==0]!=0)) return(0)

  # If it doesn't vary and matches, ignore it.
  d <- d[x.n!=0]
  x <- x[,x.n!=0,drop=FALSE]

  if(!is.null(y)){
    y <- y[,x.n!=0,drop=FALSE]
    # y, if it's given, is the constrained sample, so it's OK if it
    # doesn't vary. (E.g, the extreme case --- completely observed
    # network --- is just one configuration of statistics.)
    y.n <- effectiveSize(y)
    y.n[y.n==0] <- 1 # The actual number is irrelevant, since the cov. mat. will be 0.
  }
  x.n <- x.n[x.n!=0]
  
  v <- t(cov(x)/sqrt(x.n))/sqrt(x.n)
  if(!is.null(y)) v <- v + t(cov(y)/sqrt(y.n))/sqrt(y.n)
  chi2 <- t(d)%*%robust.inverse(v)%*%d
  pchisq(chi2,ncol(x),lower.tail=FALSE)
}

#################
# Block of code A removed from above:

##  Check for degeneracy if new network has fewer than 49999 edges
#    if(z$nedges >= 50000-1 || ergm.checkdegeneracy(statsmatrix, statsmatrix.obs, verbose=verbose)){
#      if(iteration <= control$MCMLE.maxit){
#        cat(paste("The MCMC sampler is producing degenerate samples.\n",
#                  "Try starting the algorithm at an alternative model\n",
#                  "(That is, changing the 'init' argument).\n",
#                  "I am trying something simple...\n",
#                  "The current init is:\n"))
#        print(mcmc.init)
#        mcmc.init <- 0.9*mcmc.init
#        next
#      }else{
#        cat(paste("The MCMC sampler is producing degenerate samples.\n",
#                  "Try starting the algorithm at an alternative model\n",
#                  "(That is, changing the 'init' argument).\n",
#                  "The current init is:\n"))
#        print(mcmc.init)
#        v$coef <- mcmc.init
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
#      print(summary(model$formula, basis=nw.returned)-model$target.stats)
#      if(network.naedgecount(nw) > 0){
#        cat("Summary of simulation, relative to missing network:\n")
#        a = apply(statsmatrix.miss,2,summary.statsmatrix.ergm)[4,]
#        b = sweep(apply(statsmatrix,2,summary.statsmatrix.ergm),2,a,"-")
#        print(b,scipen=6)
#        degreedist(nw.miss.returned)
#        cat("Meanstats of simulation, relative to missing network:\n")
#        print(summary(model$formula, basis=nw.miss.returned)-model$target.stats)
#        nw.returned <- network.copy(nw.miss.returned)
#      }
#    }


