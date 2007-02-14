ergm.robmon.dyn <- function(theta0, nw, model, Clist, BD, 
                        burnin, interval, proposaltype,
                        verbose=FALSE, 
                        algorithm.control=list() ){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm
  
  # The 'contol' list could have any of the following elements:
  #     phase1_n: Sample size for phase 1
  #     initial_gain:  Initial value of a for phase 2
  #     nsubphases:  Number of subphases for phase 2
  #     niterations: Initial number of iterations for first subphase of phase 2
  #     phase3_n:  Sample size for phase 3
  
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for theta0)
  n1 <- algorithm.control$phase1_n
  if(is.null(n1)) {n1 <- 7 + 3 * Clist$nparam} #default value
  eta0 <- ergm.eta(theta0, model$etamap)
  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
  MCMCparams <- list(samplesize=n1, burnin=burnin, interval=interval,
                     orig.obs=Clist$obs, meanstats=Clist$meanstats,
                     proportionbreak=algorithm.control$proportionbreak,
                     dyninterval=algorithm.control$dyninterval
                    )
  MHproposal <- list(package=algorithm.control$proposalpackage, type=proposaltype)
  cat(paste("Phase 1: ",n1,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  z <- ergm.getMCMCDynsample(nw, model, MHproposal, eta0, MCMCparams, verbose, BD)
  ubar <- apply(z$statsmatrix, 2, mean)
  Ddiag <- apply(z$statsmatrix^2, 2, mean) - ubar^2
  # This is equivalent to, but more efficient than,
  # Ddiag <- diag(t(z$statsmatrix) %*% z$statsmatrix / phase1_n - outer(ubar,ubar))
  cat("Phase 1 complete; estimated variances are:\n")
  print(Ddiag)
  
  #phase 2:  Main phase
  a <- algorithm.control$initial_gain
  if(is.null(a)) {a <- 0.1} #default value
  n_sub <- algorithm.control$nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- algorithm.control$niterations
  if(is.null(n_iter)) {n_iter <- 7 + Clist$nparam} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + Clist$nparam and a maximum of 207+Clist$nparam, with the actual 
  # number determined by the autocorrelation in the samples.
  # Thus, our default value assumes independence (for now!)
  theta <- theta0
  Ddiaginv<-1/Ddiag
  oldthetas <- NULL 
  MCMCparams$samplesize <- 1 # With samplesize=1, interval is irrelevant and burnin is crucial.
  for(subphase in 1:n_sub) {
    thetamatrix <- NULL # Will hold matrix of all theta values for this subphase
    cat(paste("Phase 2, subphase",subphase,": a=",a,",",n_iter,"iterations"))
    cat(paste(" (burnin=",MCMCparams$burnin,")\n",sep=""))
    for(iteration in 1:n_iter) {
      eta <- ergm.eta(theta, model$etamap)
      z <- ergm.getMCMCDynsample(nw, model, MHproposal, eta, MCMCparams, verbose=FALSE, BD)
      # MCMCparams$burnin should perhaps be increased here, since
      # each iteration begins from the observed network, which must be 
      # "forgotten".
      thetamatrix <- rbind(thetamatrix,theta)
      theta <- theta - a * Ddiaginv * z$statsmatrix
    }
    a <- a/2
    n_iter <- round(n_iter*2.52) # 2.52 is approx. 2^(4/3)
    thetamatrix <- rbind(thetamatrix,theta)
    theta <- apply(thetamatrix, 2, mean)
  }
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- algorithm.control$phase3_n
  if(is.null(n3)) {n3 <- 100} #default
  MCMCparams$samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
 #cat(paste(" theta=",theta,")\n",sep=""))
  eta <- ergm.eta(theta, model$etamap)
  z <- ergm.getMCMCDynsample(nw, model, MHproposal, eta, MCMCparams, verbose, BD)
# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model,
                   statsmatrix=z$statsmatrix,
                   nr.maxit=algorithm.control$nr.maxit, 
                   calc.mcmc.se=algorithm.control$calc.mcmc.se,
                   hessian=algorithm.control$hessian,
                   method=algorithm.control$method,
                   metric=algorithm.control$metric,
                   compress=algorithm.control$compress, verbose=verbose)

  endrun <- burnin+interval*(ve$samplesize-1)
  attr(ve$sample, "mcpar") <- c(burnin+1, endrun, interval)
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
                 bounddeg=BD, formula=model$formula, 
                 interval=interval, burnin=burnin, 
                 network=nw, proposaltype=proposaltype, rm.coef=theta)),
             class="ergm")
}
