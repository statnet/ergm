#  File R/vcov.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' @describeIn ergm extracts the variance-covariance matrix of
#'   parameter estimates.
#'
#' @param object {an `ergm` object.}
#' @param sources For the `vcov` method, specify whether to return
#'   the covariance matrix from the ERGM model, the estimation
#'   process, or both combined.
#'
#' @examples
#' \donttest{
#' # Extract parameter estimates as a numeric vector:
#' coef(gest)
#' # Sources of variation in parameter estimates:
#' vcov(gest, sources="model")
#' vcov(gest, sources="estimation")
#' vcov(gest, sources="all") # the default
#' }
#' @import stats
#' @export
vcov.ergm <- function(object, sources=c("all","model","estimation"), ...){
  sources <- match.arg(sources)

  src.mod <- sources %in% c("all", "model")
  src.est <- sources %in% c("all", "estimation")

  p <- nparam(object)
  
  if(src.mod){
    if(is.null(object$hessian) && is.null(object$covar)){
      object$covar <- matrix(NA, p, p)
    }
    v.mod <- NVL(object$covar, sginv(-object$hessian, tol=.Machine$double.eps^(3/4)))
    v.mod %[.|.]% ((diag(.) < 0) %|% TRUE | is.infinite(coef(object))) <- NA
    v.mod %[.|.]% object$offset <- 0
    rowcolnames(v.mod) <- param_names(object)
  }

  if(src.est){
    v.est<- NVL(object$est.cov, matrix(0, p, p))
    v.est %[.|.]% (diag(.) < 0) <- NA
    v.est %[.|.]% object$offset <- 0
    rowcolnames(v.est) <- param_names(object)
  }

  switch(sources,
         all = v.mod + v.est,
         model = v.mod,
         estimation = v.est)  
}
