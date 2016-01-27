#  File R/ergm.robmon.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
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
#   init         : the initial theta values
#   nw             : the network
#   model          : the model, as returned by <ergm.getmodel>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
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

ergm.robmon <- function(init, nw, model,
                        MHproposal,
                        verbose=FALSE, 
                        control=control.ergm() ){
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for init)
  n1 <- control$SA.phase1_n
  if(is.null(n1)) {n1 <- 7 + 3 * model$etamap$etalength} #default value
  eta0 <- ergm.eta(init, model$etamap)
  cat("Robbins-Monro algorithm with theta_0 equal to:\n")
  print(init)
#  stats <- matrix(0,ncol=model$etamap$etalength,nrow=n1)
#  stats[1,] <- model$nw.stats - model$target.stats
## stats[,]<-  rep(model$nw.stats - model$target.stats,rep(nrow(stats),model$etamap$etalength))
## control$stats <- stats
  control$MCMC.samplesize=n1
  cat(paste("Phase 1: ",n1,"iterations"))
  cat(paste(" (interval=",control$MCMC.interval,")\n",sep=""))
  z <- ergm.getMCMCsample(nw, model, MHproposal, eta0, control, verbose)
  steplength <- control$MCMLE.steplength
  # post-processing of sample statistics:  Shift each row by the
  # matrix model$nw.stats - model$target.stats, attach column names
  statsmatrix <- sweep(z$statsmatrix, 2, model$nw.stats - model$target.stats, "+")
  colnames(statsmatrix) <- model$coef.names

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
  cat("Phase 1 complete; estimated variances are:\n")
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
  if(control$parallel>0){
   control$MCMC.samplesize <- control$MCMC.samplesize*control$parallel
  }
  for(subphase in 1:n_sub) {
    thetamatrix <- NULL # Will hold matrix of all theta values for this subphase
    cat(paste("Phase 2, subphase",subphase,": a=",a,",",n_iter,"iterations"))
    cat(paste(" (burnin=",control$MCMC.burnin,")\n",sep=""))
    for(iteration in 1:n_iter) {
#cat(paste("theta:",theta,"\n"))
      eta <- ergm.eta(theta, model$etamap)
#cat(paste("eta:",eta,"\n"))

      # control$MCMC.burnin should perhaps be increased here, since
      # each iteration begins from the observed network, which must be 
      # "forgotten".
      z <- ergm.getMCMCsample(nw, model, MHproposal, eta, control, verbose=FALSE)
      # post-processing of sample statistics:  Shift each row by the
      # matrix model$nw.stats - model$target.stats, attach column names
      statsmatrix <- sweep(z$statsmatrix, 2, model$nw.stats - model$target.stats, "+")
      colnames(statsmatrix) <- model$coef.names

      thetamatrix <- rbind(thetamatrix,theta)
      statsmean <- apply(statsmatrix,2,base::mean)
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
    theta <- apply(thetamatrix, 2, base::mean)
  }
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- control$SA.phase3_n
  if(is.null(n3)) {n3 <- 20} #default
  control$MCMC.samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",control$MCMC.interval,")\n",sep=""))
  eta <- ergm.eta(theta, model$etamap)
  control$nmatrixentries = control$MCMC.samplesize * model$etamap$etalength
  z <- ergm.getMCMCsample(nw, model, MHproposal, eta, control, verbose=FALSE)
  # post-processing of sample statistics:  Shift each row by the
  # matrix model$nw.stats - model$target.stats, attach column names
  statsmatrix <- sweep(z$statsmatrix, 2, model$nw.stats - model$target.stats, "+")
  colnames(statsmatrix) <- model$coef.names

# ubar <- apply(z$statsmatrix, 2, base::mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- ginv(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

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
