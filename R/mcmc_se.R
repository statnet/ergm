#  File R/ergm.MCMCse.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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
#'   sampling.
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

  if (metric == "IS") {
    # Calculate prediction probabilities that had been used in the last MCMLE update.
    # Estimate the Hessian.
    w <- IS_weights(xsims, etaparam)
    H <- lweighted.var(gsim, w$lw)
  } else w <- IS_weights(xsims, 0)

  # Identify canonical parameters corresponding to non-offset statistics that do not vary
  novar <- diag(H) < sqrt(.Machine$double.eps)

  #  Calculate the auto-covariance of the MCMC suff. stats.
  #  and hence the MCMC s.e.
  cov.zbar <- vcov_wmean_ar(gsims, w$ws)$v

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  if(!is.null(statsmatrices.obs)){
    if (metric == "IS") {
      w.obs <- IS_weights(xsims.obs, etaparam)
      H.obs <- lweighted.var(gsim.obs, w.obs$lw)
    } else w.obs <- IS_weights(xsims.obs, 0)

    novar.obs <- diag(H.obs)<sqrt(.Machine$double.eps)
    if(any(novar&!novar.obs)) warning("Non-varying statistics in the unconstrained sample vary in the constrained sample. This should not be happening.")

    # Handle the corner case in which the constrained statistics do not vary.
    cov.zbar.obs <- if(all(novar.obs)) matrix(0, length(novar.obs), length(novar.obs))
                    else vcov_wmean_ar(gsims.obs, w.obs$ws)$v
  }else{
    cov.zbar.obs <- cov.zbar
    cov.zbar.obs[,] <- 0
    H.obs <- H
    H.obs[,] <- 0
  }

  H <- H[!novar, !novar, drop=FALSE]
  H.obs <- H.obs[!novar, !novar, drop=FALSE]

  cov.zbar <- cov.zbar[!novar, !novar, drop=FALSE]
  cov.zbar.obs <- cov.zbar.obs[!novar, !novar, drop=FALSE]


  mc.cov.offset <- matrix(0, ncol=length(theta),nrow=length(theta))

  H <- H.obs - H # Bread^-1
  cov.zbar <- cov.zbar + cov.zbar.obs # Filling

  mc.cov <- matrix(NA,ncol=length(novar),nrow=length(novar))

  if(sum(!novar)==0 || ERRVL2(srcond(-H) < .Machine$double.eps, TRUE)){
    warning("Approximate Hessian matrix is singular. Standard errors due to MCMC approximation of the likelihood may be unreliable. This is likely due to insufficient MCMC sample size or highly correlated model terms.", call.=FALSE)
  }

  mc.cov %[.,.]% !novar <- sandwich_sginv(-H, cov.zbar)
  mc.cov.offset %[.,.]% !offsettheta <- mc.cov

  rowcolnames(mc.cov.offset) <- param_names(model)

  mc.cov.offset
}

#' Time-series-adjusted variance of a weighted mean with endogenous weights
#'
#' This function assumes that the weights are endogenous and uses the
#' Delta method.
#'
#' @param x a list of matrices (potentially an [`mcmc.list`]).
#' @param w a list of vectors of row weights (unnormalized).
#'
#' @return a list with elements `m` for the mean and `v` for the
#'   variance, `infl` for inflation factor due to autocorrelation,
#'   `rank` for the rank of covariance matrix, `n` for total sample size,
#'   `neff` for total effective sample size.
#'
#' @noRd
vcov_wmean_ar <- function(x, w, tol = sqrt(.Machine$double.eps)) {
  n <- sum(lengths(w))
  xw <- map2(x, w, `*`)
  xww <- map2(xw, w, cbind)
  v <- spectrum0.mvar(xww, tol = tol)
  mw <- mean(unlist(w))
  m <- colMeans.mcmc.list(xw) / mw
  g <- rbind(diag(1, ncol(x[[1]])), -m) / mw
  v1 <- as.matrix(xTAx(g, v))
  v0 <- var.mcmc.list(x)
  e <- eigen(qrssolve(v0, v1), symmetric = TRUE, only.values = TRUE)$values
  infl <- prod(e[e / max(e) >= tol])
  structure(list(m = m, v = v1 / n,
                 infl = infl, rank = attr(v, "rank") - 1L,
                 n = n, neff = n / infl))
}

#' Importance sampling weights in a convenient format
#'
#' @param xsims a list of statistic matrices
#' @param deta the canonical parameter difference
#'
#' @return A list with elements `lw` = log-weights (vector) and `ws` =
#'   linear weights in list. The weights are *not* normalized.
#'
#' @noRd
IS_weights <- function(xsims, deta) {
  deta %[.]% is.nan <- 0 # NaN can only be offsets.
  if (all(deta == 0)) {
    sizes <- map_int(xsims, nrow)
    list(lw = numeric(sum(sizes)),
         ws = map(sizes, function(size) rep(1, size)))
  } else {
    nz <- (deta != 0) & !const_variables(xsims)
    lws <- map(xsims, function(xsim) drop(xsim[, nz, drop = FALSE] %*% deta[nz]))
    lwm <- max(unlist(lws))
    ws <- map(lws, function(lw) exp(lw - lwm))
    list(lw = unlist(lws), ws = ws)
  }
}

const_variables <- function(x) {
  map(x, apply, 2L, range) |>
    reduce(function(r1, r2) {
      rbind(pmin(r1[1L, ], r2[1L, ]),
            pmax(r1[2L, ], r2[2L, ]))
    }) |>
    diff() |>
    (`==`)(0)
}
