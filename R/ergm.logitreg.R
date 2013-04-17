#  File R/ergm.logitreg.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 2 functions for logistic regression
#        <ergm.logitreg>
#        <ergm.logisticdeviance>
#============================================================================





#############################################################################
# The <ergm.logitreg> function maximizes the log-likelihood via logistic
# regression
# 
# --PARAMETERS--
#   x        : the design matrix
#   y        : the binary outcomes, presumably the vector of edge values
#   wt       : a vector of wieghts for each case; default=rep(1,length(y));
#              this and the two params above are returned by <ergm.pl>
#   intercept: whether an intercept should be estimated (T or F);
#              default=FALSE
#   start    : initial values for the parameters to be optimized over
#   offset   : ?? (this is passed in as the 'foffset' returned by <ergm.pl>
#   maxit    : the maximum number of iterations to use; default=200
#   ...      : additional control parameters to be passed to the native
#              R routine <optim>
#
# --RETURNED--
#   fit: the best fit, as found by the <optim> function, and as a list of:
#        par         :  the best set of parameters found
#        value       :  the logistic deviance, as computed by
#                       <ergm.logisticdeviance> evaluated at 'par'
#        counts      :  a two-element vector giving the number of calls to the
#                       <ergm.locisticdeviance> function and the gradient
#                       function for the "BFGS" method; this excludes calls
#                       for computing the Hessian
#        convergence :  a convergence code; 0 indicates successful convergence,
#                       a positive number indicates various problems (see
#                       "?otpim" for details
#        message     :  additional information returned by the optimizer
#        hessian     :  a symmetric matrix giving an estimate of the Hessian
#                       evaluated at 'par'
#        cov.unscaled:  the covariance, as the robust inverse of the Hessian
#
#############################################################################

ergm.logitreg <- function(x, y, wt = rep(1, length(y)),
                          intercept = FALSE, start = rep(0, p),
                          offset=NULL, maxit=200, ...)
{
  gmin <- function(beta, X, y, w, offset) {
      eta <- (X %*% beta)+offset; p <- plogis(eta)
      -2 * matrix(w *dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% X
  }
  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  if(is.null(offset)) offset <- rep(0,length(y))
  dn <- dimnames(x)[[2]]
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x) + intercept
  if(intercept) {x <- cbind(1, x); dn <- c("(Intercept)", dn)}
  if(is.factor(y)) y <- (unclass(y) != 1)
  fit <- optim(start, ergm.logisticdeviance, gmin,
               X = x, y = y, w = wt, offset=offset,
               method = "BFGS", hessian=TRUE, control=list(maxit=maxit), ...)
  names(fit$par) <- dn
  fit$coef <- fit$par
  fit$deviance <- fit$value
  fit$iter <- fit$counts[1]
  asycov <- try(robust.inverse(fit$hessian), silent = TRUE)
  if (inherits(asycov, "try-error")) {
     asycov <- diag(1/diag(-fit$hessian))
  }
  fit$cov.unscaled <- asycov
# cat("\nCoefficients:\n"); print(fit$par)
# # R: use fit$value and fit$convergence
# cat("\nResidual Deviance:", format(fit$value), "\n")
  if(fit$convergence > 0)
      cat("\nConvergence code:", fit$convergence, "\n")
  invisible(fit)
}



ergm.logisticdeviance <- function(beta, X, y,
                            w=rep(1,length(y)), offset=rep(0,length(y))) {
      p <- plogis((X %*% beta)+offset)
      -sum(2 * w * ifelse(y, log(p), log(1-p)))
}

