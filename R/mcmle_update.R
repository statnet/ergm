#  File R/ergm.estimate.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
##################################################################################
# The <ergm.estimate> function searches for and returns a maximizer of the
# log-likelihood function. This function is observation-process capable.
#
# --PARAMETERS--
#   init          : the vector of theta parameters that produced 'statsmatrix'
#   model           : the model, as returned by <ergm_model>
#   statsmatrix     : the matrix of observed statistics that has already had the
#                     "observed statistics" vector subtracted out (i.e., the
#                     "observed stats" are assumed to be zero here)
#   statsmatrix.obs: the corresponding statsmatrix for the observation process;
#                     default=NULL
#   metric          : the name of a metric to use, as either "Likelihood",
#                     "lognormal", "Median.Likelihood", or "EF.Likelihood";
#                     default="lognormal"
#   calc.mcmc.se    : whether to calculate the standard errors induced by the
#                     MCMC algorithm; default=TRUE
#   hessainflag     : whether the Hessian matrix of the likelihood function
#                     should be computed; default=TRUE
#   verbose         : whether the progress of the estimation should be printed;
#                     default=FALSE
#   estimateonly    : whether only the estimates (vs. the estimates and the
#                     standard errors) should be calculated; default=FALSE
#
#
# --RETURNED--
#   an ergm object as a list containing several components; for details
#   see the return list in the <ergm> function header (<ergm.estimate>=^)
#
###################################################################################         

ergm.estimate<-function(init, model, statsmatrices, statsmatrices.obs=NULL,
                        metric="lognormal",
                        calc.mcmc.se=TRUE, hessianflag=TRUE,
                        verbose=FALSE,
                        metric.settings = list(),
                        steplen=1,
                        cov.type="normal",# cov.type="robust", 
                        estimateonly=FALSE, ...) {
  estimateonly <- estimateonly & !calc.mcmc.se
  # If there is an observation process to deal with, statsmatrices.obs
  # will not be NULL.
  obs <- !is.null(statsmatrices.obs)

  # Construct an offsetless map and convert init (possibly "curved"
  # parameters) to eta0 (canonical parameters)
  etamap <- model$etamap
  etamap.no <- deoffset.etamap(etamap, init)
  eta0 <- ergm.eta(init[!etamap$offsettheta], etamap.no)


  statsmatrices.orig <- statsmatrices
  statsmatrices.orig.obs <- statsmatrices.obs

  statsmean <- colMeans.mcmc.list(statsmatrices.orig)
  if(!is.null(statsmatrices.orig.obs)){
    statsmatrices.obs <- lapply.mcmc.list(statsmatrices.orig.obs, .shift_scale_points, statsmean, steplen) # I.e., shrink each point of statsmatrix.obs towards the centroid of statsmatrix.
  }else{
    statsmatrices <- sweep.mcmc.list(statsmatrices.orig, (1 - steplen) * statsmean, check.margin = FALSE)
  }
  
  statsmatrix <- as.matrix(statsmatrices)
  if(obs) statsmatrix.obs <- as.matrix(statsmatrices.obs)
  statsmatrix.orig <- as.matrix(statsmatrices.orig)
  if(obs) statsmatrix.orig.obs <- as.matrix(statsmatrices.orig.obs)

    
  # Copy and compress the stats matrices after dropping the offset terms.
  xsim <- compress_rows(as.logwmatrix(statsmatrix[,!etamap$offsetmap, drop=FALSE]))
  lrowweights(xsim) <- -log_sum_exp(lrowweights(xsim)) # I.e., divide all weights by their sum.

  xsim.orig <- compress_rows(as.logwmatrix(statsmatrix.orig[,!etamap$offsetmap, drop=FALSE]))
  lrowweights(xsim.orig) <- -log_sum_exp(lrowweights(xsim.orig)) # I.e., divide all weights by their sum.

  if(obs){
    xsim.obs <- compress_rows(as.logwmatrix(statsmatrix.obs[,!etamap$offsetmap, drop=FALSE]))
    lrowweights(xsim.obs) <- -log_sum_exp(lrowweights(xsim.obs)) 

    xsim.orig.obs <- compress_rows(as.logwmatrix(statsmatrix.orig.obs[,!etamap$offsetmap, drop=FALSE]))
    lrowweights(xsim.orig.obs) <- -log_sum_exp(lrowweights(xsim.orig.obs)) 
  }else{
    xsim.obs <- NULL
    xsim.orig.obs <- NULL
  }
  
  # It is assumed that the statsmatrix matrix has already had the
  # "observed statistics" subtracted out.  Another way to say this is
  # that when ergm.estimate is called, the "observed statistics"
  # should equal zero when measured on the scale of the statsmatrix
  # statistics.
  #' @importFrom robustbase covMcd
  if(cov.type=="robust"){
    tmp <- covMcd(decompress_rows(xsim, target.nrows=nrow(statsmatrix)))    
    av <- tmp$center
    V <- tmp$cov
  }else{
    av <- lweighted.mean(xsim, lrowweights(xsim))
    V <- lweighted.var(xsim, lrowweights(xsim))
  }
  
  xobs <- -av
  
  # Do the same recentering for the statsmatrix.obs matrix, if appropriate.
  # Note that xobs must be adjusted too.
  if(obs) {
    if(cov.type=="robust"){
      tmp <- covMcd(decompress_rows(xsim.obs, target.nrows=nrow(statsmatrix.obs)))    
      av.obs <- tmp$center
      V.obs <- tmp$cov
    }else{
      av.obs <- lweighted.mean(xsim.obs, lrowweights(xsim.obs))
      V.obs <- lweighted.var(xsim.obs, lrowweights(xsim.obs), 0)
    }

    xobs <- av.obs - av
  }
    
  # Choose appropriate loglikelihood, gradient, and Hessian functions
  # depending on metric chosen and also whether obs==TRUE
  if (verbose) { message("Using ", metric, " metric (see control.ergm function).") }

  nobs_metric <- function(name) {
    stop("Metric ", sQuote(name), " is not implemented for MLE for incompletely observed networks.")
  }

  DEP_METRIC_MAP <- list(Loglikelihood = "lognormal",
                         Median.Likelihood = "median",
                         EF.Likelihood = "naive")

  if (!is.null(new <- DEP_METRIC_MAP[[metric]])) {
    warning_once("Metric ", sQuote(metric), " has been deprecated in favor of ",
                 sQuote(new), " and may be removed in the future.", call. = FALSE)
    metric <- new
  }

  if (obs) {
    loglikelihoodfn <- switch(metric,
                              lognormal=llik.fun.obs.lognormal,
                              logtaylor = nobs_metric("logtaylor"),
                              median = llik.fun.obs.median,
                              llik.fun.obs.IS)
    gradientfn <- switch(metric,
                         lognormal = llik.grad.obs.lognormal,
                         llik.grad.obs.IS)
    Hessianfn <- switch(metric,
                        lognormal = llik.hessian.obs.lognormal,
                        llik.hessian.obs.IS)
  } else {
    loglikelihoodfn <- switch(metric,
                              lognormal=llik.fun.lognormal,
                              logtaylor=llik.fun.logtaylor,
                              median = llik.fun.median,
                              llik.fun.IS)
    gradientfn <- switch(metric,
                         lognormal = llik.grad.lognormal,
                         llik.grad.IS)
    Hessianfn <- switch(metric,
                        lognormal = llik.hessian.lognormal,
                        llik.hessian.IS)
  }

  # Get a maximal list of arguments to metrics, and filter it in
  # accordance to which ones are actually needed.
  llk_inputs <- c(list(theta = NULL,
                       xsim = xsim,
                       xsim.obs = xsim.obs,
                       eta0 = eta0,
                       etamap = etamap.no),
                  metric.settings)

  llk_inputs <- llk_inputs[
    list(loglikelihoodfn, gradientfn, Hessianfn) |> # Selected metrics
    map(formals) |>
    map(names) |>
    unlist() |>
    unique() |>
    intersect(names(llk_inputs))
  ]

  # Now find maximizer of approximate loglikelihood ratio l(eta) - l(eta0).
  # First: If we're using the lognormal approximation, the maximizer is
  # closed-form.  We can't use the closed-form maximizer if we are
  # dealing with a curved exponential family.
  if (!is.curved(model) &&
      (metric=="lognormal" || metric=="Likelihood") &&
      all(etamap$mintheta == -Inf) &&
      all(etamap$maxtheta == +Inf)) {
    if (obs) {
      if (verbose) { message("Using log-normal approx with missing (no optim)") }
      # Here, setting posd.tol=0 ensures that the matrix is
      # nonnegative-definite: it is possible for some simulated
      # statistics not to change, but it is not possible for the
      # constrained sample to have a higher variance than the
      # unconstrained.
      Lout <- list(hessian = -as.matrix(snearPD(V-V.obs,posd.tol=0)$mat))
    } else {
      if (verbose) { message("Using log-normal approx (no optim)") }
      Lout <- list(hessian = -V)
    }
    Lout$argument <- try(eta0 + ssolve(-Lout$hessian, xobs), silent = TRUE)
    # If there's an error, first try a robust matrix inverse.  This can often
    # happen if the matrix of simulated statistics does not ever change for one
    # or more statistics.
    if (inherits(Lout$argument,"try-error"))
      Lout$argument <- try(eta0 + sginv(-Lout$hessian) %*% xobs, silent = TRUE)
    # If there's still an error, use the Matrix package to try to find an 
    # alternative Hessian approximant that has no zero eigenvalues.
    if (inherits(Lout$argument, "try-error")) {
      if (obs) {
        Lout <- list(hessian = -(as.matrix(snearPD(V-V.obs)$mat)))
      }else{
        Lout <- list(hessian = -(as.matrix(snearPD(V)$mat)))
      }
      Lout$argument <- eta0 + ssolve(-Lout$hessian, xobs)
    }
    Lout$convergence <- 0 # maybe add some error-checking here to get other codes
    Lout$value <- 0.5 * crossprod(xobs, Lout$argument - eta0)
    hessianflag <- TRUE # to make sure we don't recompute the Hessian later on
  } else {
    # "guess" will be the starting point for the optim search algorithm.
    guess <- init[!etamap$offsettheta]
    mintheta <- etamap$mintheta[!etamap$offsettheta]
    maxtheta <- etamap$maxtheta[!etamap$offsettheta]
    
    trustfn <- function(theta) {
      # Check for box constraint violation.
      if (any(is.na(theta) | theta < mintheta | theta > maxtheta))
        return(list(value = -Inf))

      llk_inputs[[1L]] <- theta
      value <- do.call(loglikelihoodfn, llk_inputs)
      grad <- do.call(gradientfn, llk_inputs)
      hess <- do.call(Hessianfn, llk_inputs)
      hess[upper.tri(hess)]<-t(hess)[upper.tri(hess)]
#      message_print(value)
#      message_print(grad)
#      message_print(hess)
      list(value=value,gradient=as.vector(grad),hessian=hess)
    }

    if (verbose) message("Optimizing loglikelihood")
    #' @importFrom trust trust
    Lout <- trust(objfun = trustfn, parinit = guess, rinit = 1, rmax = 100,
                  parscale = rep(1, length(guess)), minimize = FALSE)
  }

  llk_inputs[[1L]] <- Lout$argument

  theta <- init
  theta[!etamap$offsettheta] <- Lout$argument
  names(theta) <- names(init)
  if (estimateonly) {
    # Output results as ergm-class object
    list(coefficients = setNames(theta, names(init)),
         MCMCtheta = init,
         samplesize = nrow(statsmatrix),
         loglikelihood = Lout$value,
         failure = FALSE)
  } else {
    gradient <- rep(NA, length=length(init))
    gradient[!etamap$offsettheta] <- do.call(gradientfn, llk_inputs)
    #
    #  Calculate the auto-covariance of the MCMC suff. stats.
    #  and hence the MCMC s.e.
    #
    mc.se <- rep(NA, length=length(theta))
    mc.cov <- matrix(NA, length(theta), length(theta))
    covar <- NA
    if(!hessianflag || steplen != 1)
      Lout$hessian <- do.call(Hessianfn, llk_inputs)

    invHessian <- try(sginv(-Lout$hessian, tol=.Machine$double.eps^(3/4)), TRUE)
    if(inherits(invHessian, "try-error")) invHessian <- Lout$hessian[] <- NA # Hessian non-SNND.

    covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
    covar %[.,.]% !etamap$offsettheta <- invHessian
    rowcolnames(covar) <- names(theta)
    He <- matrix(NA, ncol=length(theta), nrow=length(theta))
    He %[.,.]% !etamap$offsettheta <- Lout$hessian
    rowcolnames(He) <- names(theta)
    Lout$hessian <- He
    
    if(calc.mcmc.se){
      if (verbose) message("Starting MCMC s.e. computation.")
      mcse.metric <-
        if ((metric == "lognormal" || metric == "Likelihood") &&
            length(etamap$curved) == 0) "lognormal"
        else "IS"
      mc.cov <- ERRVL2(ergm.MCMCse(model = model, theta = theta, init = init,
                                   statsmatrices = statsmatrices,
                                   statsmatrices.obs = statsmatrices.obs,
                                   H = V, H.obs = V.obs,
                                   metric = mcse.metric),
                       matrix(NA, length(theta), length(theta), dimnames=list(names(theta),names(theta))))
    }

    # Output results as ergm-class object
    list(
      coefficients = setNames(theta, names(init)),
      sample = statsmatrices, sample.obs = statsmatrices.obs,
      iterations = Lout$counts[1], MCMCtheta = init,
      loglikelihood = Lout$value, gradient = gradient, hessian = Lout$hessian,
      covar = covar, failure = FALSE, mc.cov = mc.cov
    )
  }
}
