#  File R/ergm.logitreg.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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
                          offset=NULL, m=NULL, maxit=200, ...)
{
  if(is.null(m)){
    etamap <- identity
    etagrad <- function(theta) diag(1,length(theta),length(theta))
  }else{
    etamap <- function(theta) ergm.eta(theta,m$etamap)
    etagrad <- function(theta) ergm.etagrad(theta, m$etamap)
  }
  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  if(is.null(offset)) offset <- rep(0,length(y))
  dn <- if(is.null(m)) dimnames(x)[[2]] else coef.names.model(m, FALSE)
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x) + intercept
  if(intercept) {x <- cbind(1, x); dn <- c("(Intercept)", dn)}
  if(is.factor(y)) y <- (unclass(y) != 1)
  start[is.na(start)]<-0


  loglikelihoodfn.trust <-
    if(is.null(m)){
      function(theta, X, y, w, offset, etamap, etagrad){
        eta <- as.vector(.multiply.with.inf(X,etamap(theta))+offset)
        p <- plogis(eta)
        list(value = sum(w * ifelse(y, log(p), log1p(-p))),
             gradient = as.vector(matrix(w *dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% X %*% t(etagrad(theta))),
             hessian = -etagrad(theta) %*% crossprod(X*w*p*(1-p), X) %*% t(etagrad(theta)))
      }
    }else{
      function(theta.no, X, y, w, offset, etamap, etagrad){
        theta <- start
        theta[!m$etamap$offsettheta] <- theta.no
        
        eta <- as.vector(.multiply.with.inf(X,etamap(theta))+offset)
        p <- plogis(eta)
        list(
          value = sum(w * ifelse(y, log(p), log1p(-p))),
          gradient = as.vector((matrix(w *dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% X %*% t(etagrad(theta)))[,!m$etamap$offsettheta]),
          hessian = -(etagrad(theta) %*% crossprod(X*w*p*(1-p), X) %*% t(etagrad(theta)))[!m$etamap$offsettheta,!m$etamap$offsettheta]
        )
      }
    }

  init <- if(!is.null(m)) start[!m$etamap$offsettheta] else start
  fit <- trust(objfun=loglikelihoodfn.trust, parinit=init,
               rinit=1, 
               rmax=100, 
               parscale=rep(1,length(init)), iterlim=maxit, minimize=FALSE,
               X = x, y = y, w = wt, offset=offset, etamap=etamap, etagrad=etagrad)
  fit$coef <- fit$argument
  names(fit$coef) <- dn[if(!is.null(m)) !m$etamap$offsettheta else TRUE]
  fit$deviance <- -2*fit$value
  fit$iter <- fit$iterations
  asycov <- ginv(-fit$hessian)
  fit$cov.unscaled <- asycov
  if(!fit$converged)
      message("Trust region algorithm did not converge.")
  invisible(fit)
}



ergm.logisticdeviance <- function(theta, X, y,
                            w=rep(1,length(y)), offset=rep(0,length(y)), etamap=identity, etagrad=NULL) {
      p <- plogis(.multiply.with.inf(X,etamap(theta))+offset)
      -2*sum(w * ifelse(y, log(p), log1p(-p)))
}

