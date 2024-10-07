#  File R/ergm.stocapprox.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
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
#   s       : ergm state object
#   s.obs   : ergm state object (non NULL if there is a missing data)
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
#   verbose   : whether the MCMC sampling should be verbose (T or F);
#               default=FALSE
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.stocapprox>=%)
#
###########################################################################      

ergm.stocapprox <- function(init, s, s.obs,
                            control,
                            verbose=FALSE){
  if(!is.null(s.obs)) stop("Stochastic approximation does not support missing data at this time.")

  model <- s$model

  control <- remap_algorithm_MCMC_controls(control, "SA")

  for(ctl in c("SA.phase1_n", "SA.min_iterations", "SA.max_iterations"))
    if(is.function(control[[ctl]]))
      control[[ctl]] <- control[[ctl]](q=nparam(model, offset=FALSE), p=nparam(model, offset=FALSE, canonical=TRUE), n=network.size(s))

  message("Stochastic approximation algorithm with theta_0 equal to:")
  print(init)

  theta <- init

  control$MCMC.samplesize <- control$SA.niterations   # Now the number of Phase 2 in ergm.phase2

  for(i in 1:control$SA.nsubphases){
    control$MCMC.samplesize <- trunc(control$MCMC.samplesize*2.52)+1 # 2.52 is approx. 2^(4/3)
  }
  z <- ergm.phase12(s, theta, control, verbose=verbose)

  theta <- z$theta
  names(theta) <- names(init)
  message("Stochastic Approximation estimate:")
  message_print(theta)
  
  ## phase 3:  Estimate covariance matrix for final theta
  control$MCMC.samplesize <- control$SA.phase3_n
  message(paste("Phase 3: ",control$SA.phase3_n,"iterations"),appendLF=FALSE)
  message(paste(" (interval=",control$MCMC.interval,")",sep=""))

  # Obtain MCMC sample
  z <- ergm_MCMC_sample(z$state, control, theta=theta, verbose=max(verbose-1,0))

  if(verbose){message("Calling MCMLE Optimization...")}
  if(verbose){message("Using Newton-Raphson Step ...")}

  v <- ergm.estimate(init=theta, model=model,
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

  v$sample <- z$stats
  nws.returned <- lapply(z$networks, as.network)
  v$newnetworks <- nws.returned
  v$newnetwork <- nws.returned[[1]]
  v$coef.init <- init
  v$est.cov <- v$mc.cov
  v$mc.cov <- NULL
  v$control <- control

  v$etamap <- model$etamap
  v$MCMCflag <- TRUE
  v
}
