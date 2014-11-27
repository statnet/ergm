#  File R/ergm.stocapprox.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
############################################################################
# The <ergm.stocapprox> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is based on Snijders
# (2002), J of Social Structure and Snijders and van Duijn (2002) from
# A Festscrift for Ove Frank.  Both papers are available from Tom Snijders'
# web page:         http://stat.gamma.rug.nl/snijders/publ.htm
# The other MLE styles are found in functions <ergm.robmon>, <ergm.stepping>
# and <ergm.mainfitloop>
#
# --PARAMETERS--
#   init    : the initial theta values
#   nw        : the network
#   model     : the model, as returned by <ergm.getmodel>
#   Clist     : a list of several network and model parameters,
#               as returned by <ergm.Cprepare>
#   initialfit: an ergm object, as the initial fit
#   control: a list of parameters for controlling the MCMC sampling;
#               recognized components include
#                  'phase1_n'      'phase3_n'    'epsilon'
#                  'initial_gain'  'nsubphases'  'niterations'
#                  'nr.maxit'      'nr.reltol'   'calc.mcmc.se'
#                  'hessian'       'method'      'metric'
#                  'compress'      'trustregion' 'burnin'
#                  'interval'
#               the use of these variables is explained in the
#               <control.ergm> function header
#   MHproposal: an MHproposal object for 'nw', as returned by
#               <getMHproposal>
#   verbose   : whether the MCMC sampling should be verbose (T or F);
#               default=FALSE
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.stocapprox>=%)
#
###########################################################################      

ergm.stocapprox <- function(init, nw, model, Clist,
                            control, MHproposal,
                            verbose=FALSE){
    
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for init)
  n1 <- control$SA.phase1_n
  if(is.null(n1)) {n1 <- max(200,7 + 3 * model$etamap$etalength)} #default value
  eta0 <- ergm.eta(init, model$etamap)
  cat("Stochastic approximation algorithm with theta_0 equal to:\n")
  print(init)
  control <- c(control, list(phase1=n1,
                  stats=summary.statistics.network(model$formula, basis=nw)-model$target.stats,
                  target.stats=model$target.stats)
                 )
# cat(paste("Phase 1: ",n1,"iterations"))
# cat(paste(" (interval=",control$MCMC.interval,")\n",sep=""))
  nw.orig <- nw
  #phase 2:  Main phase
  a <- control$SA.initial_gain
  if(is.null(a)) {a <- 0.1} #default value
  n_sub <- control$SA.nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- control$SA.niterations
  if(is.null(n_iter)) {n_iter <- 7 + model$etamap$etalength} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + model$etamap$etalength and a maximum of 207+model$etamap$etalength, with the actual 
  # number determined by the autocorrelation in the samples.
  # Thus, our default value assumes independence (for now!)
  theta <- init
# Ddiaginv<-1/Ddiag
  control$MCMC.samplesize <- n_iter # Now the number of Phase 2 in ergm.phase2
  control$nsub <- n_iter # Now the number of Phase 2 sub-phases start point
  control$gain <- a # Now the number of Phase 2 sub-phases start point
# if(control$parallel>0){
#  control$MCMC.samplesize <- control$MCMC.samplesize*control$parallel
# }
  eta <- ergm.eta(theta, model$etamap)
  for(i in 1:n_sub){
    control$MCMC.samplesize <- trunc(control$MCMC.samplesize*2.52)+1 # 2.52 is approx. 2^(4/3)
  }
# cat(paste("Phase 2: a=",a,"Total Samplesize",control$MCMC.samplesize,"\n"))
# aDdiaginv <- a * Ddiaginv
  z <- ergm.phase12(nw, model, MHproposal, 
                    eta, control, verbose=TRUE)
  nw <- z$newnetwork
# toggle.dyads(nw, head = z$changed[,2], tail = z$changed[,3])
# control$maxchanges <- z$maxchanges
  theta <- z$eta
  names(theta) <- names(init)
  cat(paste(" (eta[",seq(along=theta),"] = ",paste(theta),")\n",sep=""))
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- control$SA.phase3_n
  if(is.null(n3)) {n3 <- 1000} #default
  control$MCMC.samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",control$MCMC.interval,")\n",sep=""))
#cat(paste(" (samplesize=",control$MCMC.samplesize,")\n",sep=""))
#cat(paste(" theta=",theta,")\n",sep=""))
  eta <- ergm.eta(theta, model$etamap)
#cat(paste(" (samplesize=",control$MCMC.samplesize,")\n",sep=""))
#cat(paste(" eta=",eta,")\n",sep=""))

  # Obtain MCMC sample
  z <- ergm.getMCMCsample(nw, model, MHproposal, eta0, control, verbose)
  
  # post-processing of sample statistics:  Shift each row,
  # attach column names
  statshift <- summary.statistics.network(model$formula, basis=nw) - model$target.stats
  statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
  colnames(statsmatrix) <- model$coef.names
  #v$sample <- statsmatrix
# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(init=theta, model=model,
                   statsmatrix=statsmatrix,
                   statsmatrix.obs=NULL,
                   epsilon=control$epsilon, 
                   nr.maxit=control$MCMLE.NR.maxit, 
                   nr.reltol=control$MCMLE.NR.reltol,
                   calc.mcmc.se=control$MCMC.addto.se,
                   hessianflag=control$main.hessian,
                   method=control$MCMLE.method,
                   metric=control$MCMLE.metric,
                   trustregion=control$SA.trustregion,
                   compress=control$MCMC.compress, verbose=verbose)
#
# Important: Keep R-M (pre-NR) theta
# ve$coef <- theta
#
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
                 interval=control$MCMC.interval, burnin=control$MCMC.burnin, 
                 network=nw.orig, est.cov=ve$mc.cov)),
             class="ergm")
}
