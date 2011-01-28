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
#   theta0    : the initial theta values
#   nw        : the network
#   model     : the model, as returned by <ergm.getmodel>
#   Clist     : a list of several network and model parameters,
#               as returned by <ergm.Cprepare>
#   initialfit: an ergm object, as the initial fit
#   MCMCparams: a list of parameters for controlling the MCMC sampling;
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

ergm.stocapprox <- function(theta0, nw, model, Clist,
                            MCMCparams, MHproposal,
                            verbose=FALSE){
    
  #phase 1:  Estimate diagonal elements of D matrix (covariance matrix for theta0)
  n1 <- MCMCparams$phase1_n
  if(is.null(n1)) {n1 <- max(200,7 + 3 * Clist$nstats)} #default value
  eta0 <- ergm.eta(theta0, model$etamap)
  cat("Stochastic approximation algorithm with theta_0 equal to:\n")
  print(theta0)
  names(Clist$obs) <- names(theta0)
  if(is.null(Clist$meanstats)){Clist$meanstats <- Clist$obs}
  MCMCparams <- c(MCMCparams, list(phase1=n1,
                  stats=summary.statistics.network(model$formula, basis=nw)-Clist$meanstats,
                  meanstats=Clist$meanstats)
                 )
# cat(paste("Phase 1: ",n1,"iterations"))
# cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
  nw.orig <- nw
  #phase 2:  Main phase
  a <- MCMCparams$initial_gain
  if(is.null(a)) {a <- 0.1} #default value
  n_sub <- MCMCparams$nsubphases
  if(is.null(n_sub)) {n_sub <- 4} #default value
  n_iter <- MCMCparams$niterations
  if(is.null(n_iter)) {n_iter <- 7 + Clist$nstats} #default value
  # This default value is very simplistic; Snijders would use a minimum of
  # 7 + Clist$nstats and a maximum of 207+Clist$nstats, with the actual 
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
  eta <- ergm.eta(theta, model$etamap)
  for(i in 1:n_sub){
    MCMCparams$samplesize <- trunc(MCMCparams$samplesize*2.52)+1 # 2.52 is approx. 2^(4/3)
  }
# cat(paste("Phase 2: a=",a,"Total Samplesize",MCMCparams$samplesize,"\n"))
# aDdiaginv <- a * Ddiaginv
  z <- ergm.phase12(nw, model, MHproposal, 
                    eta, MCMCparams, verbose=TRUE)
  nw <- z$newnetwork
# toggle.dyads(nw, head = z$changed[,2], tail = z$changed[,3])
# MCMCparams$maxchanges <- z$maxchanges
  theta <- z$eta
  names(theta) <- names(theta0)
  cat(paste(" (eta[",seq(along=theta),"] = ",paste(theta),")\n",sep=""))
  
  #phase 3:  Estimate covariance matrix for final theta
  n3 <- MCMCparams$phase3_n
  if(is.null(n3)) {n3 <- 1000} #default
  MCMCparams$samplesize <- n3
  cat(paste("Phase 3: ",n3,"iterations"))
  cat(paste(" (interval=",MCMCparams$interval,")\n",sep=""))
#cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
#cat(paste(" theta=",theta,")\n",sep=""))
  eta <- ergm.eta(theta, model$etamap)
#cat(paste(" (samplesize=",MCMCparams$samplesize,")\n",sep=""))
#cat(paste(" eta=",eta,")\n",sep=""))

  # Obtain MCMC sample
  MCMCparams$nmatrixentries <- Clist$nstats * MCMCparams$samplesize #Unnecessary?
  z <- ergm.getMCMCsample(Clist, MHproposal, eta0, MCMCparams, verbose)
  
  # post-processing of sample statistics:  Shift each row,
  # attach column names
  statshift <- summary.statistics.network(model$formula, basis=nw) - Clist$meanstats
  statsmatrix <- sweep(z$statsmatrix, 2, statshift, "+")
  colnames(statsmatrix) <- model$coef.names
  v$sample <- statsmatrix
# ubar <- apply(z$statsmatrix, 2, mean)
# hessian <- (t(z$statsmatrix) %*% z$statsmatrix)/n3 - outer(ubar,ubar)
# covar <- robust.inverse(covar)
  
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}

  ve<-ergm.estimate(theta0=theta, model=model,
                   statsmatrix=statsmatrix,
                   statsmatrix.miss=NULL,
                   epsilon=MCMCparams$epsilon, 
                   nr.maxit=MCMCparams$nr.maxit, 
                   nr.reltol=MCMCparams$nr.reltol,
                   calc.mcmc.se=MCMCparams$calc.mcmc.se,
                   hessianflag=MCMCparams$hessian,
                   method=MCMCparams$method,
                   metric=MCMCparams$metric,
                   trustregion=MCMCparams$trustregion,
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
                 interval=MCMCparams$interval, burnin=MCMCparams$burnin, 
                 network=nw.orig)),
             class="ergm")
}
