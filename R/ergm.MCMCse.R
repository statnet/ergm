#  File R/ergm.MCMCse.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#############################################################################
# The <ergm.MCMCse> function computes and returns the MCMC standard errors
#
# --PARAMETERS--
#   theta           :  the vector of theta coefficients
#   init          :  the vector of initial theta coefficients
#   statsmatrix     :  the matrix of network statistics
#   statsmatrix.obs :  the matrix of network statistics on the constrained network
#   model           :  the model, as returned by <ergm_model>
#
# --RETURNED--
#   the variance of the MCMC sampling as a list containing:
#     mc.se : the vector of MCMC standard error estimates for each theta parameter
#     mc.cov: the MCMC covariance matrix of the theta parameters
#
################################################################################

ergm.MCMCse<-function(theta, init, statsmatrix, statsmatrix.obs,
                      model) {
  # Not sure why this is necessary, but:
  names(theta) <- names(init)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(init, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap
  etaparam <- (eta-eta0)[!offsetmap]

  # Center statmatrix (and statsmatrix.obs, if applicable)
  av <- apply(statsmatrix, 2, mean)
# av <- apply(statsmatrix,2,median)
  xsim <- sweep(statsmatrix, 2, av, "-")
  gsim <- ergm.estfun(xsim, theta, model)
  xobs <- -av
  if(!is.null(statsmatrix.obs)){
   av.obs <- apply(statsmatrix.obs, 2, mean)
#  av.obs <- apply(statsmatrix.obs, 2, median)
   xsim.obs <- sweep(statsmatrix.obs, 2, av.obs,"-")
   gsim.obs <- ergm.estfun(xsim.obs, theta, model)
   xsim.obs <- xsim.obs[,!offsetmap, drop=FALSE]
   xobs <- av.obs-av
  }
  xobs <- xobs[!offsetmap]
  xsim <- xsim[,!offsetmap, drop=FALSE]

  # Calculate prediction probabilities that had been used in the last MCMLE update.
  basepred <- xsim %*% etaparam
  prob <- max(basepred)
  prob <- exp(basepred - prob)
  prob <- prob/sum(prob)

  # Estimate the Hessian. Do as much as possible on the scale of the estimating functions.
  E <- apply(sweep(gsim, 1, prob, "*"), 2, sum)
  htmp <- sweep(sweep(gsim, 2, E, "-"), 1, sqrt(prob), "*")
  H <- crossprod(htmp, htmp)

  #  Calculate the auto-covariance of the MCMC suff. stats.
  #  and hence the MCMC s.e.
  cov.zbar <- spectrum0.mvar(gsim) * sum(prob^2)
  imp.factor <- sum(prob^2)*length(prob)

  # Identify canonical parameters corresponding to non-offset statistics that do not vary
  novar <- rep(TRUE, nrow(H))
  novar <- diag(H) < sqrt(.Machine$double.eps)

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  if(!is.null(statsmatrix.obs)){
    obspred <- xsim.obs %*% etaparam
    prob.obs <- exp(obspred - max(obspred))
    prob.obs <- prob.obs/sum(prob.obs)
    E.obs <- apply(sweep(gsim.obs, 1, prob.obs, "*"), 2, sum)
    htmp.obs <- sweep(sweep(gsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
    H.obs <- crossprod(htmp.obs, htmp.obs)

    cov.zbar.obs <- spectrum0.mvar(gsim.obs) * sum(prob.obs^2)
    imp.factor.obs <- sum(prob.obs^2)*length(prob.obs)

    novar.obs <- diag(H.obs)<sqrt(.Machine$double.eps)
    if(any(novar&!novar.obs)) warning("Non-varying statistics in the unconstrained sample vary in the constrained sample. This should not be happening.")
  }else{
    cov.zbar.obs <- cov.zbar
    cov.zbar.obs[,] <- 0
    H.obs <- H
    H.obs[,] <- 0
    imp.factor.obs <- NULL
  }

  H <- H[!novar, !novar, drop=FALSE]
  H.obs <- H.obs[!novar, !novar, drop=FALSE]

  cov.zbar <- cov.zbar[!novar, !novar, drop=FALSE]
  cov.zbar.obs <- cov.zbar.obs[!novar, !novar, drop=FALSE]


  mc.cov.offset <- matrix(0, ncol=length(theta),nrow=length(theta))

  H <- H.obs - H # Bread^-1
  cov.zbar <- cov.zbar + cov.zbar.obs # Filling

  mc.cov <- matrix(NA,ncol=length(novar),nrow=length(novar))

  if(sum(!novar)==0 || inherits(try(solve(H,tol=1e-20),silent=TRUE),"try-error")){
    warning("Approximate Hessian matrix is singular. Standard errors due to MCMC approximation of the likelihood cannot be evaluated. This is likely due to insufficient MCMC sample size or highly correlated model terms.", call.=FALSE)
  }else{
    mc.cov0 <- solve(H, cov.zbar, tol=1e-20)
    mc.cov0 <- solve(H, t(mc.cov0), tol=1e-20)
    mc.cov[!novar,!novar] <- mc.cov0
  }

  mc.cov.offset[!offsettheta,!offsettheta] <- mc.cov

  colnames(mc.cov.offset) <- names(theta)
  rownames(mc.cov.offset) <- names(theta)

  attr(mc.cov.offset, "imp.factor") <- imp.factor
  attr(mc.cov.offset, "imp.factor.obs") <- imp.factor.obs
  
  mc.cov.offset
}
