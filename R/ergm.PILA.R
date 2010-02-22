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
