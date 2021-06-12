#  File R/ergm.MCMCse.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' Compute the MCMC standard errors
#'
#' @param theta the vector of theta coefficients
#' @param init the vector of initial theta coefficients
#' @param statsmatrices an [`mcmc.list`] of  network statistics
#' @param statsmatrices.obs an [`mcmc.list`] of constrained network statistics
#' @param H the Hessian matrix (lognormal metric only)
#' @param H.obs the Hessian matrix on the constrained network (lognormal metric only)
#' @param model the [`ergm_model`]
#'
#' @return The variance matrix of the parameter estimates due to MCMC
#'   sampling, with attributes `"imp.factor"` and `"imp.factor.obs"`
#'   giving the additional variation due to importance sampling; the
#'   latter are always 1 for lognormal metric.
#' @noRd
ergm.MCMCse <- function(model, theta, init, statsmatrices, statsmatrices.obs,
                        H, H.obs, metric = c("IS", "lognormal")) {
  metric <- match.arg(metric)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(init, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap
  if(metric == "IS") etaparam <- (eta-eta0)[!offsetmap]

  # Center statmatrix (and statsmatrices.obs, if applicable)
  av <- colMeans.mcmc.list(statsmatrices)
# av <- apply(statsmatrices,2,median)
  xsims <- sweep.mcmc.list(statsmatrices, av, "-")
  gsims <- ergm.estfun(xsims, theta=theta, model=model)
  xobs <- -av
  xsims <- xsims[,!offsetmap, drop=FALSE]
  xsim <- as.matrix(xsims)
  gsim <- as.matrix(gsims)

  if(!is.null(statsmatrices.obs)){
   av.obs <- colMeans.mcmc.list(statsmatrices.obs)
#  av.obs <- apply(statsmatrices.obs, 2, median)
   xsims.obs <- sweep.mcmc.list(statsmatrices.obs, av.obs,"-")
   gsims.obs <- ergm.estfun(xsims.obs, theta=theta, model=model)
   xsims.obs <- xsims.obs[,!offsetmap, drop=FALSE]
   xsim.obs <- as.matrix(xsims.obs)
   gsim.obs <- as.matrix(gsims.obs)

   xobs <- av.obs-av
  }
  xobs <- xobs[!offsetmap]

  if(metric == "IS") {
    # Calculate prediction probabilities that had been used in the last MCMLE update.
    basepred <- xsim %*% etaparam
    prob <- max(basepred)
    prob <- exp(basepred - prob)
    prob <- prob/sum(prob)

    # Estimate the Hessian. Do as much as possible on the scale of the estimating functions.
    E <- apply(sweep(gsim, 1, prob, "*"), 2, sum)
    htmp <- sweep(sweep(gsim, 2, E, "-"), 1, sqrt(prob), "*")
    H <- crossprod(htmp, htmp)
  } else prob <- rep.int(1 / nrow(xsim), nrow(xsim))

  #  Calculate the auto-covariance of the MCMC suff. stats.
  #  and hence the MCMC s.e.
  cov.zbar <- spectrum0.mvar(gsims) * sum(prob^2)
  imp.factor <- sum(prob^2)*length(prob)

  # Identify canonical parameters corresponding to non-offset statistics that do not vary
  novar <- rep(TRUE, nrow(H))
  novar <- diag(H) < sqrt(.Machine$double.eps)

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  if(!is.null(statsmatrices.obs)){
    if(metric == "IS") {
      obspred <- xsim.obs %*% etaparam
      prob.obs <- exp(obspred - max(obspred))
      prob.obs <- prob.obs/sum(prob.obs)
      E.obs <- apply(sweep(gsim.obs, 1, prob.obs, "*"), 2, sum)
      htmp.obs <- sweep(sweep(gsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
      H.obs <- crossprod(htmp.obs, htmp.obs)
    } else prob.obs <- rep.int(1 / nrow(xsim.obs), nrow(xsim.obs))

    cov.zbar.obs <- spectrum0.mvar(gsims.obs) * sum(prob.obs^2)
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

  rownames(mc.cov.offset) <- colnames(mc.cov.offset) <- param_names(model)

  attr(mc.cov.offset, "imp.factor") <- imp.factor
  attr(mc.cov.offset, "imp.factor.obs") <- imp.factor.obs
  
  mc.cov.offset
}
