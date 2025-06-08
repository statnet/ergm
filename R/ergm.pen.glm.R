#  File R/ergm.pen.glm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#===================================================================
# This file contains the following 3 files for penalized glm fits
#             <ergm.pen.glm>
#             <model.matrix.pen.glm>
#===================================================================

.BernLogLik <- function(p, y, w=1){
  stopifnot(all(y%in%c(0,1)))
  sum(ifelse(y, log(p), log1p(-p))*w)
}



###############################################################################
# The <ergm.pen.glm> function calculates and returns an ergm fit using a
# a penalized glm approach
#
# --PARAMETERS--
#   formula: a formula 'y ~ x + ..', presumably y as the 'zy' returned by 
#            <ergm.pl> and x as 'xmat' returned by <ergm.pl>
#   data   : the dataframe or environment in which to find the variables
#            contained in 'formula'; defaults to the formula's environment.
#   alpha  : ??, this code merely returns 'alpha', and no caller provides
#            an alpha; default=0.05
#   maxit  : the maximum number of iterations to use in this fitting;
#            default=25
#   maxhs  : the maximum number of half steps to use; default=5
#   epsilon: the value used to return the fit; when the sum of changes ?
#            falls below this, the fit is returned; if 'maxit' iterations
#            are realized, this is ignored; default=1e-4
#   maxstep: ??; default=10
#   start  : an optional starting vector of theta coefficients;
#            default=NULL
#   weights: the weights corresponding to y and x of 'formula'; default=NULL
#
#
# --RETURNED--
#   fit: a 'pen.glm' object as a list containing:
#     coefficients     : the theta coefficients maximizing the penalized
#                        log-likelihood
#     alpha            : the 'alpha' above 
#     var              : the covariance matrix
#     df               : the degrees of freedom
#     loglik           : the log-likelihood corresponding to 'coefficients'
#     iter             : the number of iterations used, as either 'maxit' or
#                        the number needed to achieve convergence according
#                        to 'epsilon'
#     terms            : the names of the coefficients
#     formula          : the 'formula' above
#     data             : the 'data' above
#     model.matrix     : the design matrix
#     method           : "pen.glm.fit"
#     linear.predictors: the predicted y as 'model.matrix'*'coefficients'
#
##############################################################################

ergm.pen.glm <- function(formula,
  data, alpha = 0.05,
  maxit = 25, maxhs = 5, epsilon = 0.0001, maxstep = 10, 
  start=NULL,
  weights=NULL)
{
  if(missing(data)) data <- environment(formula)
  y <- as.vector(model.extract(model.frame(formula, data = data), "response"))
  n <- length(y)
  x <- model.matrix(formula, data = data)  # Model-Matrix

#  aaa <- rep(seq(along=weights),weights)
#  x <- as.matrix(x[aaa,])
#  y <- y[aaa]
#  n <- length(y)
#  weights <- y-y+1

  k <- ncol(x)  # Anzahl Effekte

  int <- 0
  coltotest <-1:k

  if(missing(weights)){weights <- rep(1,length=n)}
  beta <- c(log((sum(y*weights)/sum((1-y)*weights))),
            rep(0, k - 1))
  if(!missing(start) && !is.null(start) && ncol(x)==length(start) && !is.na(start)){
    beta[1] <- beta[1] - sum((x %*% start)*weights)
  }
  iter <- 0
  pi <- expit(drop(x %*% beta))
  loglik <- .BernLogLik(pi, y, weights)
  XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*")  #### X' (W ^ 1/2)
  Fisher <- crossprod(XW2)  #### X' W  X
  loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
  repeat {
   iter <- iter + 1
   XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
   Fisher <- crossprod(XW2)  #### X' W  X
   covs <- sginv(Fisher, tol=.Machine$double.eps^(3/4))  ### (X' W  X) ^ -1
#  H <- crossprod(XW2, covs) %*% XW2
#  H <- XW2 %*% covs %*% t(XW2)
   diagH <- pi
   for(i in seq(along=diagH)){
    diagH[i] <- xTAx(XW2[i,], covs)
   }
#  U.star <- crossprod(x, y - pi)
#  U.star <- crossprod(x, (y - pi)*weights)
#  U.star <- crossprod(x, (y - pi + diagH * (0.5 - pi))*weights)
   U.star <- crossprod(x, (y - pi)*weights + diagH * (0.5 - pi))
   delta <- as.vector(covs %*% U.star)
   mx <- max(abs(delta))/maxstep
   if(mx > 1){
    delta <- delta/mx
   }
   beta <- beta + delta
   loglik.old <- loglik
   for(halfs in 1:maxhs) {
## Half-Steps
    pi <- expit(drop(x %*% beta))
    loglik <- .BernLogLik(pi, y, weights)
    XW2 <- sweep(x, 1, (weights*(pi * (1 - pi)))^0.5, "*") #### X' (W ^ 1/2)
    Fisher <- crossprod(XW2) #### X' W  X
    loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    if(loglik > loglik.old) break
    beta <- beta - delta * 2^( - halfs) ##beta-A enderung verkleinern
   }
   if(iter == maxit | sum(abs(delta)) <= epsilon){break}
  }

  fit <- list(coefficients = beta, alpha = alpha, var = covs, 
              df = (k-int), 
              loglik = loglik,
              deviance = -2*loglik,
              iter = iter, n = n, 
              terms = colnames(x), y = y, 
              formula = as.formula(formula), call=match.call(),
              data = data, 
              model.matrix = x, 
              method = "pen.glm.fit",
              linear.predictors=as.vector(x %*% beta))
  class(fit) <- c("pen.glm")
# vars <- diag(covs)
  fit
}

#' @noRd
model.matrix.pen.glm <- function(object, ...)
{
	object$model.matrix
}
