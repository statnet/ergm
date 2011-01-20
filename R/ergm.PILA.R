######################################################################
# The <ergm.PILA> function estimates an ergm, but uniquely does so by
# using the PILA sampler to create the initial theta vector, before
# the usual estimation process continues
#
# --PARAMETERS--
#   theta0    : the initial theta values used to draw a PILA sample
#               and subsequently find the usual 'theta0'
#   nw        : the network
#   model     : the model, as returned by <ergm.getmodel>
#   Clist     : a list of several network and model parameters,
#               as returned by <ergm.Cprepare>
#   MHproposal: an MHproposal object for 'nw', as returned by
#               <getMHproposal>
#   MCMCparams: a list of parameters for controlling the MCMC sampling;
#               recognized components include:
#     samplesize     : the size of the sample to collect
#     interval       : the number of proposals to disregard between
#                      samples
#     burnin         : the number of proposals to initially disregard
#                      for the burn-in period
#     PILA.steplength: ??
#     PILA.gamma     : ??
#     parallel       : the number of threads in which to run the sampling
#     packagenames   : names of packages; this is only relevant if
#                      "ergm" is given
#     Clist.dt       : this is a Clist, similar to that returned by
#                      <ergm.Cprepare>, but this is for fitting dynamic
#                      models
#     Clist.miss     : a corresponding 'Clist' for the network of missing 
#                      edges, as returned by <ergm.design>
#   control   : a list of control parameters for estimation as
#               returned by <control.ergm>
#   verbose   : whether the MCMC sampling should be verbose (T or F);
#               default=FALSE
#   
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.PILA>= +);
#
######################################################################

ergm.PILA <- function(theta0, nw, model, Clist,
                      MHproposal,MCMCparams, control,
                      verbose=FALSE){
  if(verbose)cat("PILA algorithm with theta_0 equal to:\n")
  print(theta0)
  stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)
  stats[1,] <- Clist$obs - Clist$meanstats

  MCMCparams$stats<-stats
  
  if(verbose)cat("Running PILA:\n")

  z <- ergm.runPILAsampler(nw, model, MHproposal, theta0, MCMCparams, verbose)

  theta<-z$etamatrix[dim(z$etamatrix)[1],]
  if(verbose)cat("Theta estimated:",theta,"\n")

  hist<-z
  z <- ergm.getMCMCsample.parallel(nw, model, MHproposal, theta, MCMCparams, verbose)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model,
                   statsmatrix=z$statsmatrix,
                   statsmatrix.miss=NULL,
                   nr.maxit=control$nr.maxit, 
                   nr.reltol=control$nr.reltol,
                   calc.mcmc.se=control$calc.mcmc.se,
                   hessianflag=control$hessian,
                   method=control$method,
                   metric=control$metric,
                   compress=control$compress, verbose=verbose)

  endrun <- with(MCMCparams,burnin+interval*(ve$samplesize-1))
  attr(ve$sample, "mcpar") <- with(MCMCparams,c(burnin+1, endrun, interval))
  attr(ve$sample, "class") <- "mcmc"
  ve$null.deviance <- 2*network.dyadcount(nw)*log(2)
  ve$mle.lik <- -ve$null.deviance/2 + ve$loglikelihood
# The next is the right one to uncomment
# ve$mcmcloglik <- ve$mcmcloglik - network.dyadcount(nw)*log(2)

  # From ergm.estimate:
  #    structure(list(coef=theta, sample=statsmatrix, 
                      # iterations=iteration, mcmcloglik=mcmcloglik,
                      # MCMCtheta=theta0, 
                      # loglikelihood=loglikelihood, gradient=gradient,
                      # covar=covar, samplesize=samplesize, failure=FALSE,
                      # mc.se=mc.se, acf=mcmcacf,
                      # fullsample=statsmatrix.all),
                  # class="ergm") 
  structure(c(ve, list(newnetwork=nw, 
                 theta.original=theta0,
                 PILA.coef=theta,
                       PILA.hist=hist,
                 interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                 network=nw)),
             class="ergm")
}
