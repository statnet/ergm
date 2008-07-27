ergm.robmon.dyn <- function(theta0, nw, model.form, model.diss, Clist, 
                            gamma0,
                            MCMCparams, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm

  eta0 <- ergm.eta(theta0, model.form$etamap)

  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
  eta0 <- ergm.eta(theta0, model.form$etamap)

  z <- ergm.phase12.dyn(nw, Clist$meanstats, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta0, gamma0, MCMCparams, verbose=verbose)

  eta <- z$eta
  MCMCparams$samplesize<-MCMCparams$RobMon.phase3n
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta, gamma0,
                             MCMCparams, verbose)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=eta, model=model.form,
                    statsmatrix=z$statsmatrix.form,
                    statsmatrix.miss=NULL,
                    nr.maxit=MCMCparams$nr.maxit, 
                    nr.reltol=MCMCparams$nr.reltol, 
                    trustregion=MCMCparams$trustregion, 
                    calc.mcmc.se=MCMCparams$calc.mcmc.se,
                    hessian=MCMCparams$hessian,
                    method=MCMCparams$method,
                    metric=MCMCparams$metric,
                    compress=MCMCparams$compress, verbose=verbose)
  #
  # Important: Keep R-M (pre-NR) theta
  ve$coef <- eta
  #
  endrun <- MCMCparams$burnin+MCMCparams$interval*(ve$samplesize-1)
  attr(ve$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
  attr(ve$sample, "class") <- "mcmc"
  ve$null.deviance <- 2*network.dyadcount(nw)*log(2)
  ve$mle.lik <- -ve$null.deviance/2 + ve$loglikelihood
  ve$mcmcloglik <- ve$mcmcloglik - network.dyadcount(nw)*log(2)

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
                       interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                       network=nw)),
            class="ergm")
}
