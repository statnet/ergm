ergm.robmon.dyn <- function(theta0, nw, model.form, model.diss, Clist, BD, 
                        gamma0,
                        MCMCparams, MHproposal.form, MHproposal.diss,
                        verbose=FALSE){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm
  
  # The 'algorithm.contol' list could have any of the following elements:
  #     phase1_n: Sample size for phase 1
  #     initial_gain:  Initial value of a for phase 2
  #     nsubphases:  Number of subphases for phase 2
  #     niterations: Initial number of iterations for first subphase of phase 2
  #     phase3_n:  Sample size for phase 3
  
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for theta0)
  n1 <- MCMCparams$phase1_n
  if(is.null(n1)) {n1 <- max(200,7 + 3 * Clist$nparam)} #default value
  eta0 <- ergm.eta(theta0, model.form$etamap)
  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
  names(Clist$obs) <- names(theta0)
  if(is.null(Clist$meanstats)){Clist$meanstats <- Clist$obs}
  stats <- matrix(0,ncol=Clist$nparam,nrow=MCMCparams$samplesize)
  stats[1,] <- summary.statistics.network(model.form$formula, basis=nw)-Clist$meanstats
  MCMCparams <- c(MCMCparams, list(samplesize=100, phase1=n1,
                  meanstats=Clist$meanstats,
                  stats=stats)
                 )
# cat(paste("Phase 1: ",n1,"iterations"))
# cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  nw.orig <- nw
# z <- ergm.getMCMCDynsample(nw, model.form, model.diss, MHproposal, 
#                            eta0, MCMCparams, verbose, BD)
# toggle.dyads(nw, tail = z$changed[,2], head = z$changed[,3])
# nw <- network.update(nw, z$newedgelist)
# MCMCparams$orig.obs <- summary(model.form$formula)
# MCMCparams$maxchanges <- z$maxchanges
# ubar <- apply(z$statsmatrix, 2, mean)
# Ddiag <- apply(z$statsmatrix^2, 2, mean) - ubar^2
  # This is equivalent to, but more efficient than,
  # Ddiag <- diag(t(z$statsmatrix) %*% z$statsmatrix / phase1_n - outer(ubar,ubar))
# cat("Phase 1 complete; estimated variances are:\n")
# print(Ddiag)
# print(1/Ddiag)
  
  #phase 2:  Main phase
  a <- MCMCparams$initial_gain
  if(is.null(a)) {a <- 0.1} #default value
  n_sub <- MCMCparams$nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- MCMCparams$niterations
  if(is.null(n_iter)) {n_iter <- 7 + Clist$nparam} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + Clist$nparam and a maximum of 207+Clist$nparam, with the actual 
  # number determined by the autocorrelation in the samples.
  # Thus, our default value assumes independence (for now!)
  theta <- theta0
# Ddiaginv<-1/Ddiag
  MCMCparams$samplesize <- n_iter # Now the number of Phase 2 in ergm.phase2
  MCMCparams$nsub <- n_iter # Now the number of Phase 2 sub-phases start point
  MCMCparams$gain <- a # Now the number of Phase 2 sub-phases start point
# if(MCMCparams$parallel>0){
#  MCMCparams$samplesize <- MCMCparams$samplesize*MCMCparams$parallel
# }
  eta <- ergm.eta(theta, model.form$etamap)
  for(i in 1:n_sub){
    MCMCparams$samplesize <- trunc(MCMCparams$samplesize*2.52)+1 # 2.52 is approx. 2^(4/3)
  }
# cat(paste("Phase 2: a=",a,"Total Samplesize",MCMCparams$samplesize,"\n"))
# aDdiaginv <- a * Ddiaginv
  z <- ergm.phase12.dyn(nw, model.form, model.diss, MHproposal.form, MHproposal.diss,
                        eta, gamma0, MCMCparams, verbose=TRUE, BD)
  cat("\n Number changed: ");cat(paste(nrow(z$changed)));cat("\n Edges:")
# cat(paste(summary(nw ~ edges)));cat("\n")
# nw <- network.update(nw, z$newedgelist)
# toggle.dyads(nw, tail = z$changed[,2], head = z$changed[,3])
  stats[1,] <- summary.statistics.network(model.form$formula, basis=nw)-Clist$meanstats
  MCMCparams$stats <- stats
# print( MCMCparams$orig.obs)
# MCMCparams$maxchanges <- z$maxchanges
  theta <- z$eta
  names(theta) <- names(theta0)
  cat(paste(" (eta= ",paste(theta),")\n",sep=""))
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- MCMCparams$phase3_n
  if(is.null(n3)) {n3 <- 1000} #default
  MCMCparams$samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
#cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
#cat(paste(" theta=",theta,")\n",sep=""))
  eta <- ergm.eta(theta, model.form$etamap)
#cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
#cat(paste(" eta=",eta,")\n",sep=""))
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss, 
                             MHproposal.form, MHproposal.diss, eta, gamma0,
                             MCMCparams, verbose, BD)
  MCMCparams$maxchanges <- z$maxchanges
# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model.form,
                   statsmatrix=z$statsmatrix,
                   statsmatrix.miss=NULL,
                   nr.maxit=MCMCparams$nr.maxit, 
                   trustregion=MCMCparams$trustregion, 
                   calc.mcmc.se=MCMCparams$calc.mcmc.se,
                   hessian=MCMCparams$hessian,
                   method=MCMCparams$method,
                   metric=MCMCparams$metric,
                   compress=MCMCparams$compress, verbose=verbose)
#
# Important: Keep R-M (pre-NR) theta
# ve$coef <- theta
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
                 bounddeg=BD, formula=model.form$formula, 
                 interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                 network=nw.orig, proposalname=MHproposal$name)),
             class="ergm")
}

ergm.robmon.dyn.orig <- function(theta0, nw, model.form, model.diss, Clist, BD, 
                        gamma0,
                        burnin, interval, MHproposal,
                        verbose=FALSE, 
                        algorithm.control=list() ){
  # This is based on Snijders (2002), J of Social Structure
  # and Snijders and van Duijn (2002) from A Festscrift for Ove Frank
  # Both papers are available from Tom Snijders' web page: 
  #          http://stat.gamma.rug.nl/snijders/publ.htm
  
  # The 'algorithm.contol' list could have any of the following elements:
  #     phase1_n: Sample size for phase 1
  #     initial_gain:  Initial value of a for phase 2
  #     nsubphases:  Number of subphases for phase 2
  #     niterations: Initial number of iterations for first subphase of phase 2
  #     phase3_n:  Sample size for phase 3
  
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for theta0)
  n1 <- algorithm.control$phase1_n
  if(is.null(n1)) {n1 <- max(200,7 + 3 * Clist$nparam)} #default value
  eta0 <- ergm.eta(theta0, model.form$etamap)
  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
  MCMCparams <- list(samplesize=n1, burnin=burnin, interval=interval,
                     orig.obs=Clist$obs, meanstats=Clist$meanstats,
                     parallel=algorithm.control$parallel,
                     dyninterval=algorithm.control$dyninterval,
                     maxchanges=algorithm.control$maxchanges
                    )
  cat(paste("Phase 1: ",n1,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss, MHproposal.form, MHproposal.diss, eta0, gamma0, MCMCparams, verbose, BD)
  MCMCparams$maxchanges <- z$maxchanges
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
  MCMCparams$samplesize <- 1 # With samplesize=1, interval is irrelevant and burnin is crucial.
  if(MCMCparams$parallel>0){
   MCMCparams$samplesize <- MCMCparams$samplesize*MCMCparams$parallel
  }
  for(subphase in 1:n_sub) {
    thetamatrix <- NULL # Will hold matrix of all theta values for this subphase
    cat(paste("Phase 2, subphase",subphase,": a=",a,",",n_iter,"iterations"))
    cat(paste(" (burnin=",MCMCparams$burnin,")\n",sep=""))
    for(iteration in 1:n_iter) {
      eta <- ergm.eta(theta, model.form$etamap)
      if(verbose){
       cat(paste(" (iteration = ",iteration," of ",n_iter,")\n",sep=""))
      }
      z <- ergm.getMCMCDynsample(nw, model.form, model.diss, MHproposal.form, MHproosal.diss, 
                                 eta, gamma0, MCMCparams, verbose=FALSE, BD)
      MCMCparams$maxchanges <- z$maxchanges
      # MCMCparams$burnin should perhaps be increased here, since
      # each iteration begins from the observed network, which must be 
      # "forgotten".
      thetamatrix <- rbind(thetamatrix,theta)
      statsmean <- apply(z$statsmatrix,2,mean)
      theta <- theta - a * Ddiaginv * statsmean
    }
    a <- a/2
    n_iter <- round(n_iter*2.52) # 2.52 is approx. 2^(4/3)
    thetamatrix <- rbind(thetamatrix,theta)
    theta <- apply(thetamatrix, 2, mean)
  }
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- algorithm.control$phase3_n
  if(is.null(n3)) {n3 <- 1000} #default
  MCMCparams$samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
 #cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
 #cat(paste(" theta=",theta,")\n",sep=""))
  eta <- ergm.eta(theta, model.form$etamap)
  z <- ergm.getMCMCDynsample(nw, model.form, model.diss, MHproposal.form, MHproposal.diss, eta, gamma0, MCMCparams, verbose, BD)
  MCMCparams$maxchanges <- z$maxchanges
# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model.form,
                   statsmatrix=z$statsmatrix,
                   statsmatrix.miss=NULL,
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
  structure(c(ve, list(newnetwork=network.update(nw, z$newedgelist), 
                 theta.original=theta0,
                 bounddeg=BD, formula=model.form$formula, 
                 interval=interval, burnin=burnin,
                 network=nw, proposalname.form=MHproposal.form$name, proposalname.diss=MHproposal.diss$name)),
             class="ergm")
}
