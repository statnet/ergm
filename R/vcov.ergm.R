#  File R/vcov.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

#' @describeIn ergm extracts estimated model coefficients.
#' 
#' @param object {an `ergm` object.}
#' @examples
#' \donttest{
#' # Extract parameter estimates as a numeric vector:
#' coef(gest)
#' }
#' @import stats
#' @importFrom stats coef
#' @export
coef.ergm <- function(object, ...){object$coef}

#' @describeIn ergm An \emph{alias} for
#'   \code{ergm}.
#' @importFrom stats coefficients
#' @export
coefficients.ergm <- coef.ergm

#' @describeIn ergm extracts the variance-covariance matrix of
#'   parameter estimates.
#' 
#' @param sources For the `vcov` method, specify whether to return
#'   the covariance matrix from the ERGM model, the estimation
#'   process, or both combined.
#'
#' @examples
#' \donttest{
#' # Sources of variation in parameter estimates:
#' vcov(gest, sources="model")
#' vcov(gest, sources="estimation")
#' vcov(gest, sources="all") # the default
#' }
#' @importFrom stats vcov
#' @export
vcov.ergm <- function(object, sources=c("all","model","estimation"), ...){
  sources <- match.arg(sources)

  src.mod <- sources %in% c("all", "model")
  src.est <- sources %in% c("all", "estimation")

  p <- length(object$coef)
  
  if(src.mod){
    if(is.null(object$hessian) && is.null(object$covar)){
      object$covar <- matrix(NA, p, p)
    }
    v.mod <- NVL(object$covar, ginv(-object$hessian))
    v.mod[is.na(diag(v.mod))|diag(v.mod)<0|is.infinite(object$coef),] <- NA
    v.mod[,is.na(diag(v.mod))|diag(v.mod)<0|is.infinite(object$coef)] <- NA
    v.mod[object$offset,] <- 0
    v.mod[,object$offset] <- 0
     colnames(v.mod) <- rownames(v.mod) <- names(object$coef)
  }

  if(src.est){
    v.est<- NVL(object$est.cov, matrix(0, p, p))
    v.est[diag(v.est)<0,] <- NA
    v.est[,diag(v.est)<0] <- NA
    v.est[object$offset,] <- 0
    v.est[,object$offset] <- 0
    colnames(v.est) <- rownames(v.est) <- names(object$coef)
  }

  switch(sources,
         all = v.mod + v.est,
         model = v.mod,
         estimation = v.est)  
}
