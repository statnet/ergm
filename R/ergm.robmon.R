#  File R/ergm.robmon.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
############################################################################
# The <ergm.robmon> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is based on Snijders
# (2002), J of Social Structure  and Snijders and van Duijn (2002) from
# A Festscrift for Ove Frank.  Both papers are available from Tom Snijders'
# web page:           https://stat.gamma.rug.nl/snijders/publ.htm
# (The other MLE styles are found in functions <ergm.stocapprox>,
# <ergm.stepping> and <ergm.mainfitloop>
#
#
# --PARAMETERS--
#   init         : the initial theta values
#   nw             : the network
#   model          : the model, as returned by <ergm_model>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
#   burnin         : the number of proposals to disregard before sampling
#                    begins
#   interval       : the number of proposals between sampled statistics;
#   initialfit     : an ergm object, as the initial fit
#   proposal     : an proposal object for 'nw', as returned by
#                    <getproposal>
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

ergm.robmon <- function(init, nw, model,
                        proposal,
                        verbose=FALSE, 
                        control=control.ergm() ){
  # Start cluster if required (just in case we haven't already).
  ergm.getCluster(control, max(verbose-1,0))

  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for init)
  n1 <- control$SA.phase1_n
  if(is.null(n1)) {n1 <- 7 + 3 * model$etamap$etalength} #default value
  eta0 <- ergm.eta(init, model$etamap)
  message("Robbins-Monro algorithm with theta_0 equal to:")
  print(init)
#  stats <- matrix(0,ncol=model$etamap$etalength,nrow=n1)
#  stats[1,] <- model$nw.stats - model$target.stats
## stats[,]<-  rep(model$nw.stats - model$target.stats,rep(nrow(stats),model$etamap$etalength))
## control$stats <- stats
  control$MCMC.samplesize=n1
  message(paste("Phase 1: ",n1,"iterations"),appendLF=FALSE)
  message(paste(" (interval=",control$MCMC.interval,")",sep=""))
  z <- ergm_MCMC_sample(nw, model, proposal, control, eta=eta0, verbose=max(verbose-1,0))
  steplength <- control$MCMLE.steplength
  # post-processing of sample statistics:  Shift each row by the
  # matrix model$nw.stats - model$target.stats, attach column names
  statsmatrix <- sweep(as.matrix(z$stats), 2, model$nw.stats - model$target.stats, "+")
  colnames(statsmatrix) <- param_names(model,canonical=TRUE)

  if(steplength<1){
    statsmean <- apply(statsmatrix,2,base::mean)
    statsmatrix <- sweep(statsmatrix,2,(1-steplength*0.1)*statsmean,"-")
  }
# ubar <- apply(z$statsmatrix, 2, base::mean)
# Ddiag <- apply(z$statsmatrix^2, 2, base::mean) - ubar^2
# Ddiag <- apply(z$statsmatrix, 2, stats::var)
  Ddiag <- apply(statsmatrix^2, 2, base::mean)
  # This is equivalent to, but more efficient than,
  # Ddiag <- diag(t(z$statsmatrix) %*% z$statsmatrix / phase1_n - outer(ubar,ubar))
  message("Phase 1 complete; estimated variances are:")
  print(Ddiag)
#browser()
# require(covRobust)
# Ddiag <- diag(cov.nnve(z$statsmatrix))
  #phase 2:  Main phase
  a <- control$SA.initial_gain
  if(is.null(a)) {a <- 0.1/control$MCMLE.steplength} #default value
  n_sub <- control$SA.nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- control$SA.niterations
  if(is.null(n_iter)) {n_iter <- 7 + model$etamap$etalength} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + model$etamap$etalength and a maximum of 207+model$etamap$etalength, with the actual 
  # number determined by the autocorrelation in the samples.
  # Thus, our default value assumes independence (for now!)
  theta <- init
  oldthetas <- NULL 
  control$MCMC.samplesize <- 10 # With samplesize=1, interval is irrelevant and burnin is crucial.
  

  control$MCMC.samplesize <- control$MCMC.samplesize*nthreads(control)

  for(subphase in 1:n_sub) {
    thetamatrix <- NULL # Will hold matrix of all theta values for this subphase
    message(paste("Phase 2, subphase",subphase,": a=",a,",",n_iter,"iterations"), appendLF=FALSE)
    message(paste(" (burnin=",control$MCMC.burnin,")",sep=""))
    for(iteration in 1:n_iter) {
#message(paste("theta:",theta,""))
      eta <- ergm.eta(theta, model$etamap)
#message(paste("eta:",eta,""))

      # control$MCMC.burnin should perhaps be increased here, since
      # each iteration begins from the observed network, which must be 
      # "forgotten".
      z <- ergm_MCMC_sample(nw, model, proposal, control, eta=eta, verbose=max(verbose-1,0))
      # post-processing of sample statistics:  Shift each row by the
      # matrix model$nw.stats - model$target.stats, attach column names
      statsmatrix <- sweep(as.matrix(z$stats), 2, model$nw.stats - model$target.stats, "+")
      colnames(statsmatrix) <- param_names(model,canonical=TRUE)

      thetamatrix <- rbind(thetamatrix,theta)
      statsmean <- apply(statsmatrix,2,base::mean)
      if(steplength<1 && subphase < n_sub ){
        statsmean <- steplength*statsmean
      }
      
      Ddiaginv<-1/Ddiag
#message(paste("Ddiaginv:",Ddiaginv,""))
#     theta <- theta - a * Ddiaginv * statsmatrix
      theta <- theta - a * Ddiaginv * statsmean
    }
message(paste("theta new:",theta,""))
    a <- a/2
    n_iter <- round(n_iter*2.52) # 2.52 is approx. 2^(4/3)
    thetamatrix <- rbind(thetamatrix,theta)
    theta <- apply(thetamatrix, 2, base::mean)
  }
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- control$SA.phase3_n
  if(is.null(n3)) {n3 <- 20} #default
  control$MCMC.samplesize <- n3
  message(paste("Phase 3: ",n3,"iterations"),appendLF=FALSE)
  message(paste(" (interval=",control$MCMC.interval,")",sep=""))
  eta <- ergm.eta(theta, model$etamap)
  control$nmatrixentries = control$MCMC.samplesize * model$etamap$etalength
  z <- ergm_MCMC_sample(nw, model, proposal, control, eta=eta, verbose=max(verbose-1,0))
  # post-processing of sample statistics:  Shift each row by the
  # matrix model$nw.stats - model$target.stats, attach column names
  statsmatrix <- sweep(as.matrix(z$stats), 2, model$nw.stats - model$target.stats, "+")
  colnames(statsmatrix) <- param_names(model,canonical=TRUE)

# ubar <- apply(z$statsmatrix, 2, base::mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- ginv(covar)
  
  if(verbose){message("Calling MCMLE Optimization...")}
  if(verbose){message("Using Newton-Raphson Step ...")}

  ve<-ergm.estimate(init=theta, model=model,
                   statsmatrix=statsmatrix,
                   statsmatrix.obs=NULL,
                   nr.maxit=control$MCMLE.NR.maxit, 
                   nr.reltol=control$MCMLE.NR.reltol,
                   calc.mcmc.se=control$MCMC.addto.se,
                   hessianflag=control$main.hessian,
                   method=control$MCMLE.method,
                   metric=control$MCMLE.metric,
                   compress=control$MCMC.compress, verbose=verbose)

  ve$sample <- ergm.sample.tomcmc(ve$sample, control)
# The next is the right one to uncomment
# ve$mcmcloglik <- ve$mcmcloglik - network.dyadcount(nw)*log(2)

  # From ergm.estimate:
  #    structure(list(coef=theta, sample=statsmatrix, 
                      # iterations=iteration, mcmcloglik=mcmcloglik,
                      # MCMCtheta=init, 
                      # loglikelihood=loglikelihood, gradient=gradient,
                      # covar=covar, samplesize=samplesize, failure=FALSE,
                      # mc.se=mc.se, acf=mcmcacf,
                      # fullsample=statsmatrix.all),
                  # class="ergm") 
  structure(c(ve, list(newnetwork=nw, 
                 theta.original=init,
                 rm.coef=theta,
                 interval=control$MCMC.interval, burnin=control$MCMC.burnin, 
                 network=nw, est.cov=ve$mc.cov)),
             class="ergm")
}
