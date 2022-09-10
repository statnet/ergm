#  File R/ergm.stocapprox.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
############################################################################
# The <ergm.stocapprox> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is based on Snijders
# (2002), J of Social Structure and Snijders and van Duijn (2002) from
# A Festscrift for Ove Frank.  Both papers are available from Tom Snijders'
# web page:         https://stat.gamma.rug.nl/snijders/publ.htm
#
# --PARAMETERS--
#   init    : the initial theta values
#   nw        : the network
#   model     : the model, as returned by <ergm_model>
#   control: a list of parameters for controlling the MCMC sampling;
#               recognized components include
#                  'phase1_n'      'phase3_n'    'epsilon'
#                  'initial_gain'  'nsubphases'  'niterations'
#                  'nr.maxit'      'nr.reltol'   'calc.mcmc.se'
#                  'hessian'       'method'      'metric'
#                  'compress'      'burnin'
#                  'interval'
#               the use of these variables is explained in the
#               <control.ergm> function header
#   proposal: an proposal object for 'nw', as returned by
#               <getproposal>
#   verbose   : whether the MCMC sampling should be verbose (T or F);
#               default=FALSE
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.stocapprox>=%)
#
###########################################################################      

ergm.stocapprox <- function(init, nw, model,
                            control, proposal,
                            verbose=FALSE){

  control <- remap_algorithm_MCMC_controls(control, "SA")

  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for init)
  n1 <- control$SA.phase1_n
  if(is.null(n1)) {n1 <- max(200,7 + 3 * model$etamap$etalength)} #default value
  message("Stochastic approximation algorithm with theta_0 equal to:")
  print(init)
  control <- within(control, {
    phase1 <- n1
    stats <- model$nw.stats - NVL(model$target.stats,model$nw.stats)
    target.stats <- model$target.stats
  })
# message(paste("Phase 1: ",n1,"iterations"))
# message(paste(" (interval=",control$MCMC.interval,")",sep=""))
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
  for(i in 1:n_sub){
    control$MCMC.samplesize <- trunc(control$MCMC.samplesize*2.52)+1 # 2.52 is approx. 2^(4/3)
  }
# message(paste("Phase 2: a=",a,"Total Samplesize",control$MCMC.samplesize,""))
# aDdiaginv <- a * Ddiaginv
  s <- ergm_state(nw, model=model, proposal=proposal, stats = summary(model, nw) - NVL(model$target.stats,model$nw.stats))
  z <- ergm.phase12(s, 
                    theta, control, verbose=TRUE)
  nw <- z$newnetwork
# toggle.dyads(nw, head = z$changed[,2], tail = z$changed[,3])
# control$maxchanges <- z$maxchanges
  theta <- z$theta
  names(theta) <- names(init)
  message(paste(" (theta[",seq(along=theta),"] = ",paste(theta),")",sep=""))
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- control$SA.phase3_n
  if(is.null(n3)) {n3 <- 1000} #default
  control$MCMC.samplesize <- n3
  message(paste("Phase 3: ",n3,"iterations"),appendLF=FALSE)
  message(paste(" (interval=",control$MCMC.interval,")",sep=""))
#message(paste(" (samplesize=",control$MCMC.samplesize,")",sep=""))
#message(paste(" theta=",theta,")",sep=""))
  eta <- ergm.eta(theta, model$etamap)
#message(paste(" (samplesize=",control$MCMC.samplesize,")",sep=""))
#message(paste(" eta=",eta,")",sep=""))

  # Obtain MCMC sample
  z <- ergm_MCMC_sample(z$state, control, theta=theta, verbose=max(verbose-1,0))
  
#v$sample <- stats
# ubar <- apply(z$stats, 2, mean)
# hessian <- (t(z$stats) %*% z$stats)/n3 - outer(ubar,ubar)
# covar <- ginv(covar)
  
  if(verbose){message("Calling MCMLE Optimization...")}
  if(verbose){message("Using Newton-Raphson Step ...")}

  ve<-ergm.estimate(init=theta, model=model,
                   statsmatrices=mcmc.list(as.mcmc(z$stats)),
                   statsmatrices.obs=NULL,
                   epsilon=control$epsilon, 
                   nr.maxit=control$MCMLE.NR.maxit, 
                   nr.reltol=control$MCMLE.NR.reltol,
                   calc.mcmc.se=control$MCMC.addto.se,
                   hessianflag=control$main.hessian,
                   method=control$MCMLE.method,
                   metric=control$MCMLE.metric,
                   verbose=verbose)
#
# Important: Keep R-M (pre-NR) theta
# ve$coefficients <- theta
#
  ve$sample <- ergm.sample.tomcmc(ve$sample, control)
# The next is the right one to uncomment
# ve$mcmcloglik <- ve$mcmcloglik - network.dyadcount(nw)*log(2)

  # From ergm.estimate:
  #    structure(list(coefficients=theta, sample=mcmc.list(as.mcmc(stats)), 
                      # iterations=iteration, mcmcloglik=mcmcloglik,
                      # MCMCtheta=init, 
                      # loglikelihood=loglikelihood, gradient=gradient,
                      # covar=covar, samplesize=samplesize, failure=FALSE,
                      # mc.se=mc.se, acf=mcmcacf,
                      # fullsample=stats.all),
                  # class="ergm") 
  structure(c(ve, list(newnetwork=nw, 
                 theta.original=init,
                 interval=control$MCMC.interval, burnin=control$MCMC.burnin, 
                 network=nw.orig, est.cov=ve$mc.cov)),
             class="ergm")
}
