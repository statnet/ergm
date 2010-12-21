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
#     samplesize   : the number of MCMC sampled networks
#     maxit        : the maximum number of iterations to use
#     Clist.miss   : the 'Clist' for the network of missing edges, as
#                    returned by <ergm.design>
#     epsilon      : ??, this is passed to <ergm.estimate>, which ignores it;
#                    also, this is used in place of the 'epsilon' argument for
#                    function
#   MHproposal     : an MHproposal object for 'nw', as returned by
#                    <getMHproposal>
#   MHproposal.miss: an MHproposal object for the missing network of'nw',
#                    as returned by <getMHproposal>
#   verbose        : whether the MCMC sampling should be verbose (T or F);
#                    default=FALSE
#   sequential     : whether to update the network returned in
#                    'v$newnetwork'; if the network has missing edges,
#                    this is ignored; default=TRUE
#   estimate       : whether to optimize the theta0 coefficients via
#                    <ergm.estimate>; default=TRUE
#   ...            : additional parameters that may be passed from within;
#                    all are ignored
#
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.mainfitloop>=*);
#      note that if the model is degenerate, only 'coef' and 'sample' are
#      returned; if 'estimate'=FALSE, the MCMC and se variables will be
#      NA or NULL
#
###########################################################################

ergm.mainfitloop <- function(theta0, nw, model, Clist,
                             initialfit, 
                             MCMCparams, 
                             MHproposal, MHproposal.miss,
                             verbose=FALSE,
                             sequential=TRUE,
                             estimate=TRUE, ...) {
  # Calculate the amount by which all of the MCMC statistics should be adjusted
  # to account for the fact that they are all calculated relative to the
  # observed network.  Unless otherwise specified, the value of meanstats
  # is simply the observed statistics, which means statshift equals zero
  # most of the time.
  statshift <- Clist$obs - Clist$meanstats
  MCMCparams$meanstats <- Clist$meanstats
  MCMCparams$nmatrixentries = MCMCparams$samplesize * Clist$nstats

  #  while(any(mcmc.precision*asyse < mc.se, na.rm=TRUE) && iteration <= maxit){
  iteration <- 0
  finished <- FALSE
  # mcmc.theta0 will change at each iteration.  It is the value that is used
  # to generate the MCMC samples.  theta0 will never change.
  mcmc.theta0 <- theta0
  while(!finished){
	  iteration <- iteration + 1
    cat("Iteration ",iteration," of at most ", MCMCparams$maxit,
        " with parameter: \n", sep="")
    print(mcmc.theta0)

    # Obtain MCMC sample
    mcmc.eta0 <- ergm.eta(mcmc.theta0, model$etamap)
    z <- ergm.getMCMCsample(Clist, MHproposal, mcmc.eta0, MCMCparams, verbose)
    
    # post-processing of sample statistics:  Shift each row by the
    # matrix Clist$obs - Clist$meanstats, attach column names
    statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
    colnames(statsmatrix) <- model$coef.names
    newnw <- network.update(nw, z$newedgelist, "edgelist")

    ##  Does the same, if missing edges:
	  if(MCMCparams$Clist.miss$nedges > 0){
      z <- ergm.getMCMCsample(Clist, MHproposal.miss, mcmc.eta0, MCMCparams, verbose)
      statsmatrix.miss <- sweep(z$statsmatrix, 2, statshift, "+")
      colnames(statsmatrix.miss) <- model$coef.names
      if(verbose){cat("Back from constrained MCMC...\n")}
    }else{
      statsmatrix.miss <- NULL
      if(verbose){cat("Back from unconstrained MCMC...\n")}
      if (sequential) {
        #      nw.obs <- summary(model$formula, basis=newnw)
        #      namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
      }
    }

    # Removed block A of code here.  (See end of file.)
    
    if(!estimate){
      if(verbose){cat("Skipping optimization routines...\n")}
      l <- list(coef=mcmc.theta0, mc.se=rep(NA,length=length(theta0)),
                sample=statsmatrix, sample.miss=statsmatrix.miss,
                iterations=1, MCMCtheta=mcmc.theta0,
                loglikelihood=NA, #mcmcloglik=NULL, 
                mle.lik=NULL,
                gradient=rep(NA,length=length(theta0)), #acf=NULL,
                samplesize=MCMCparams$samplesize, failure=TRUE,
                newnetwork = newnw)
      return(structure (l, class="ergm"))
    }
    finished <- iteration >= MCMCparams$maxit
    # Use estimateonly=TRUE if this is not the last iteration.
    v<-ergm.estimate(theta0=mcmc.theta0, model=model,
                     statsmatrix=statsmatrix, 
                     statsmatrix.miss=statsmatrix.miss, 
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
    mcmc.theta0 <- v$coef
  } # end of main loop
  
  # FIXME:  We should not be "tacking on" extra list items to the 
  # object returned by ergm.estimate.  Instead, it is more transparent
  # if we build the output object (v) from scratch, of course using 
  # some of the info returned from ergm.estimate.
  v$sample <- statsmatrix
  v$burnin <- MCMCparams$burnin
  v$samplesize <- MCMCparams$samplesize
  v$interval <- MCMCparams$interval
  v$network <- nw
  v$newnetwork <- newnw
  v$theta.original <- theta0
  v$mplefit <- initialfit
  v$parallel <- MCMCparams$parallel

  # Removed block B of code here.  (See end of file.)

  endrun <- MCMCparams$burnin+MCMCparams$interval*(MCMCparams$samplesize-1)
  attr(v$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
  attr(v$sample, "class") <- "mcmc"
  v$null.deviance <- 2*network.dyadcount(nw)*log(2)
  v$mle.lik <- initialfit$mle.lik + abs(v$loglikelihood)
  v$etamap <- model$etamap
  v
}

#################
# Block of code A removed from above:

##  Check for degeneracy if new network has fewer than 49999 edges
#    if(NROW(z$newedgelist) >= 50000-1 || ergm.checkdegeneracy(statsmatrix, statsmatrix.miss, verbose=verbose)){
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
#    
#    if(verbose){
#      cat(paste("The density of the returned network is",
#                network.density(newnw),"\n"))
#      cat(paste("The density of the original network is",
#                network.density(nw),"\n"))
#      cat("Summary of simulation, relative to observed network:\n")
#      print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
#      degreedist(newnw)
#      cat("Meanstats of simulation, relative to observed network:\n")
#      print(summary(model$formula, basis=newnw)-Clist$meanstats)
#    }
    
#################
# Block of code B removed from above:

#  if(!v$failure & !any(is.na(v$coef))){
#    #     asyse <- sqrt(diag(robust.inverse(-v$hessian)))
#    #     asyse <- try(sqrt(diag(robust.inverse(cov(statsmatrix)))))
#    asyse <- mc.se
#    options(warn=-1)
#    #     options(warn=2)
#    if(is.null(v$covar)){
#      asyse[names(v$coef)] <- sqrt(diag(robust.inverse(-v$hessian)))
#    }else{
#      asyse[names(v$coef)] <- sqrt(diag(v$covar))
#    }
#    options(warn=0)
#  }

