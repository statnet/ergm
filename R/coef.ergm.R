#  File R/coef.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' Extract Model Fit Coefficients and Uncertainty Estimates
#' 
#' \code{coef} extracts model coefficients from objects returned by
#' the \code{\link{ergm}} function.
#' 
#' @param object {an object for which the extraction of model coefficients is
#'     meaningful.}
#' @param ... {other arguments.}
#' 
#' @return Coefficients extracted from the model object \code{object}.
#' 
#' @seealso \code{\link{fitted.values}} and \code{\link{residuals}} for related methods;
#'   \code{\link{glm}}, \code{\link{lm}} for model fitting.
#' 
#' @examples
#' data(florentine)
#' fit <- ergm(flomarriage ~ edges + concurrent)
#' coef(fit)
#' 
#' @keywords regression models
#' @export
coef.ergm <- function(object, ...){object$coef}

#' @rdname coef.ergm
#'
#' @description
#' \code{coefficients} is an \emph{alias} for \code{ergm}.
#' @export
coefficients.ergm <- coef.ergm

#' @rdname coef.ergm
#'
#' @description
#' \code{vcov} extracts the variance-covariance matrix of parameter
#'   estimates.
#' 
#' @param sources {Specify whether to return the covariance matrix
#'   from the ERGM model, the estimation process, or both combined.}
#'
#' @examples
#' vcov(fit, sources="model")
#' vcov(fit, sources="estimation")
#' vcov(fit, sources="all") # the default
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
