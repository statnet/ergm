############################################################################
# The <ergm.robmon> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is based on Snijders
# (2002), J of Social Structure  and Snijders and van Duijn (2002) from
# A Festscrift for Ove Frank.  Both papers are available from Tom Snijders'
# web page:           http://stat.gamma.rug.nl/snijders/publ.htm
# (The other MLE styles are found in functions <ergm.stocapprox>,
# <ergm.stepping> and <ergm.mainfitloop>
#
#
# --PARAMETERS--
#   theta0         : the initial theta values
#   nw             : the network
#   model          : the model, as returned by <ergm.getmodel>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
#   Clist.miss     : a 'Clist' for the network of missing edges, as
#                    returned by <ergm.design>
#   burnin         : the number of proposals to disregard before sampling
#                    begins
#   interval       : the number of proposals between sampled statistics;
#   initialfit     : an ergm object, as the initial fit
#   MHproposal     : an MHproposal object for 'nw', as returned by
#                    <getMHproposal>
#   verbose        : whether the MCMC sampling should be verbose (T or F);
#                    default=FALSE
#   control        : a list of parameters for controlling the fitting
#                    process, as returned by <control.ergm>; in
#                    particular, the following components are recognized:
#                     'phase1_n'      'parallel'    'steplength'
#                     'initial_gain'  'nsubphases'  'niterations'
#                     'nr.maxit'      'nr.reltol'   'calc.mcmc.se'
#                     'hessian'       'method'      'metric'
#                     'compress'
#
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.robmon>=&);
#
###########################################################################      

ergm.robmon <- function(theta0, nw, model, Clist, Clist.miss,
                        burnin, interval, MHproposal,
                        verbose=FALSE, 
                        control=control.ergm() ){
  
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for theta0)
  n1 <- control$phase1_n
  if(is.null(n1)) {n1 <- 7 + 3 * Clist$nstats} #default value
  eta0 <- ergm.eta(theta0, model$etamap)
  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(theta0)
#  stats <- matrix(0,ncol=Clist$nstats,nrow=n1)
#  stats[1,] <- Clist$obs - Clist$meanstats
## stats[,]<-  rep(Clist$obs - Clist$meanstats,rep(nrow(stats),Clist$nstats))
## MCMCparams$stats <- stats
  MCMCparams <- c(control, 
                  list(samplesize=n1, burnin=burnin, interval=interval,
                       nmatrixentries = n1* Clist$nstats, #stats=stats, 
                       parallel=control$parallel, Clist.miss=Clist.miss))
  cat(paste("Phase 1: ",n1,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  z <- ergm.getMCMCsample(Clist, MHproposal, eta0, MCMCparams, verbose)
  steplength <- control$steplength
  # post-processing of sample statistics:  Shift each row by the
  # matrix Clist$obs - Clist$meanstats, attach column names
  statsmatrix <- sweep(z$statsmatrix, 2, Clist$obs - Clist$meanstats, "+")
  colnames(statsmatrix) <- model$coef.names

  if(steplength<1){
    statsmean <- apply(statsmatrix,2,mean)
    statsmatrix <- sweep(statsmatrix,2,(1-steplength*0.1)*statsmean,"-")
  }
# ubar <- apply(z$statsmatrix, 2, mean)
# Ddiag <- apply(z$statsmatrix^2, 2, mean) - ubar^2
# Ddiag <- apply(z$statsmatrix, 2, var)
  Ddiag <- apply(statsmatrix^2, 2, mean)
  # This is equivalent to, but more efficient than,
  # Ddiag <- diag(t(z$statsmatrix) %*% z$statsmatrix / phase1_n - outer(ubar,ubar))
  cat("Phase 1 complete; estimated variances are:\n")
  print(Ddiag)
#browser()
# require(covRobust)
# Ddiag <- diag(cov.nnve(z$statsmatrix))
  #phase 2:  Main phase
  a <- control$initial_gain
  if(is.null(a)) {a <- 0.1/control$steplength} #default value
  n_sub <- control$nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- control$niterations
  if(is.null(n_iter)) {n_iter <- 7 + Clist$nstats} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + Clist$nstats and a maximum of 207+Clist$nstats, with the actual 
  # number determined by the autocorrelation in the samples.
  # Thus, our default value assumes independence (for now!)
  theta <- theta0
  oldthetas <- NULL 
  MCMCparams$samplesize <- 10 # With samplesize=1, interval is irrelevant and burnin is crucial.
  if(MCMCparams$parallel>0){
   MCMCparams$samplesize <- MCMCparams$samplesize*MCMCparams$parallel
  }
  for(subphase in 1:n_sub) {
    thetamatrix <- NULL # Will hold matrix of all theta values for this subphase
    cat(paste("Phase 2, subphase",subphase,": a=",a,",",n_iter,"iterations"))
    cat(paste(" (burnin=",MCMCparams$burnin,")\n",sep=""))
    for(iteration in 1:n_iter) {
#cat(paste("theta:",theta,"\n"))
      eta <- ergm.eta(theta, model$etamap)
#cat(paste("eta:",eta,"\n"))

      # MCMCparams$burnin should perhaps be increased here, since
      # each iteration begins from the observed network, which must be 
      # "forgotten".
      MCMCparams$nmatrixentries = MCMCparams$samplesize * Clist$nstats
      z <- ergm.getMCMCsample(Clist, MHproposal, eta, MCMCparams, verbose=FALSE)
      # post-processing of sample statistics:  Shift each row by the
      # matrix Clist$obs - Clist$meanstats, attach column names
      statsmatrix <- sweep(z$statsmatrix, 2, Clist$obs - Clist$meanstats, "+")
      colnames(statsmatrix) <- model$coef.names

      thetamatrix <- rbind(thetamatrix,theta)
      statsmean <- apply(statsmatrix,2,mean)
      if(steplength<1 && subphase < n_sub ){
        statsmean <- steplength*statsmean
      }
      
      Ddiaginv<-1/Ddiag
#cat(paste("Ddiaginv:",Ddiaginv,"\n"))
#     theta <- theta - a * Ddiaginv * statsmatrix
      theta <- theta - a * Ddiaginv * statsmean
    }
cat(paste("theta new:",theta,"\n"))
    a <- a/2
    n_iter <- round(n_iter*2.52) # 2.52 is approx. 2^(4/3)
    thetamatrix <- rbind(thetamatrix,theta)
    theta <- apply(thetamatrix, 2, mean)
  }
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- control$phase3_n
  if(is.null(n3)) {n3 <- 20} #default
  MCMCparams$samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  eta <- ergm.eta(theta, model$etamap)
  MCMCparams$nmatrixentries = MCMCparams$samplesize * Clist$nstats
  z <- ergm.getMCMCsample(Clist, MHproposal, eta, MCMCparams, verbose=FALSE)
  # post-processing of sample statistics:  Shift each row by the
  # matrix Clist$obs - Clist$meanstats, attach column names
  statsmatrix <- sweep(z$statsmatrix, 2, Clist$obs - Clist$meanstats, "+")
  colnames(statsmatrix) <- model$coef.names

# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model,
                   statsmatrix=statsmatrix,
                   statsmatrix.miss=NULL,
                   nr.maxit=control$nr.maxit, 
                   nr.reltol=control$nr.reltol,
                   calc.mcmc.se=control$calc.mcmc.se,
                   hessianflag=control$hessian,
                   method=control$method,
                   metric=control$metric,
                   compress=control$compress, verbose=verbose)

  endrun <- burnin+interval*(ve$samplesize-1)
  attr(ve$sample, "mcpar") <- c(burnin+1, endrun, interval)
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
                 rm.coef=theta,
                 interval=interval, burnin=burnin, 
                 network=nw)),
             class="ergm")
}
