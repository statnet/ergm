#  File R/ergm.logitreg.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#============================================================================
# This file contains the following 2 functions for logistic regression
#        <ergm.logitreg>
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
#        value       :  the logistic deviance
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
                          intercept = FALSE, start = NULL,
                          offset=NULL, m=NULL, maxit=200, verbose=FALSE, ...){

  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  if(is.null(offset)) offset <- dbl_along(y)
  if(intercept){
    x <- cbind(1, x)
    dn <- c("(Intercept)", dn)
  }
  if(is.factor(y)) y <- (unclass(y) != 1)

  if(is.null(m)){
    etamap <- identity
    etagrad <- function(theta) diag(1,length(theta),length(theta))
    dn <- dimnames(x)[[2]]
    if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
    p <- ncol(x)
  }else{
    etamap <- function(theta) ergm.eta(theta,m$etamap)
    etagrad <- function(theta) ergm.etagrad(theta, m$etamap)
    dn <- param_names(m, FALSE)
    p <- nparam(m, FALSE)
  }

  NVL(start) <- rep(NA, p)
  start %[f]% (\(x) is.na(x) & !is.nan(x)) <- rnorm(sum(is.na(start)), 0, sqrt(.Machine$double.eps))

  loglikelihoodfn.trust <-
    if(is.null(m)){
      function(theta, X, y, w, offset, etamap, etagrad){
        eta <- as.vector(mat_by_coef(X, etamap(theta)) + offset)
        Xgradt <- X %*% t(etagrad(theta))
        p <- expit(eta)
        o <- list(value = sum(w * ifelse(y, log(p), log1p(-p))),
             gradient = as.vector(matrix(w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% Xgradt),
             hessian = -crossprod(Xgradt, w*p*(1-p)*Xgradt))
        if(verbose){
          message("theta:")
          message_print(theta)
          message("result:")
          message_print(o)
        }
        o
      }
    }else{
      function(theta.no, X, y, w, offset, etamap, etagrad){
        theta <- start
        theta[!m$etamap$offsettheta] <- theta.no

        # Check for box constraint violation.
        if (any(is.na(theta.no)) ||
            any(theta.no < m$etamap$mintheta[!m$etamap$offsettheta]) ||
            any(theta.no > m$etamap$maxtheta[!m$etamap$offsettheta])) {
          o <- list(value = -Inf)
        } else {
          eta <- as.vector(mat_by_coef(X, etamap(theta)) + offset)
          ## FIXME: This can certainly be precomputed, rather than
          ## done every iteration.
          if (any(etanan <- is.nan(eta))) {
            eta <- eta[!etanan]
            X <- X[!etanan, , drop = FALSE]
            w <- w[!etanan]
            y <- y[!etanan]
          }
          Xgradt <- X %*% t(etagrad(theta))
          p <- expit(eta)
          o <- list(
            value = sum(w * ifelse(y, log(p), log1p(-p))),
            gradient = as.vector((matrix(w * dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% Xgradt)[,!m$etamap$offsettheta,drop=FALSE]),
            hessian = -crossprod(Xgradt, w*p*(1-p)*Xgradt)[!m$etamap$offsettheta,!m$etamap$offsettheta,drop=FALSE]
          )
        }
        if(verbose){
          message("theta:")
          message_print(theta)
          message("result:")
          message_print(o)
        }
        o
      }
    }

  init <- if(is.null(m)) start else .constrain_init(m, start)[!m$etamap$offsettheta]
  
  fit <- trust(objfun=loglikelihoodfn.trust, parinit=init,
               rinit=1, 
               rmax=100, 
               parscale=rep(1,length(init)), iterlim=maxit, minimize=FALSE,
               X = x, y = y, w = wt, offset=offset, etamap=etamap, etagrad=etagrad)
  fit$coefficients <- fit$argument
  names(fit$coefficients) <- dn[if(!is.null(m)) !m$etamap$offsettheta else TRUE]
  fit$deviance <- -2*fit$value
  fit$iter <- fit$iterations
  asycov <- sginv(-fit$hessian, tol=.Machine$double.eps^(3/4))
  fit$cov.unscaled <- asycov
  if(!fit$converged)
      message("Trust region algorithm did not converge.")
  fit
}
