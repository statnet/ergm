#  File R/rank_test.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#excerpted from gof.ergm.R:

##FIXME: Add support for curved ERGMs.
##FIXME: Merge with gof.ergm and/or with ergm itself.
##FIXME: Publish this somewhere.

#' A lack-of-fit test for ERGMs
#'
#' A simple test reporting the sample quantile of the observed
#' network's probability in the distribution under the MLE. This is a
#' conservative p-value for the null hypothesis of the observed
#' network being a draw from the distribution of interest.
#'
#' @param x an [ergm()] object.
#' @param plot if `TRUE`, plot the empirical distribution.
#'
#' @return The sample quantile of the observed network's probability
#'   among the predicted.
#'
#' @export
rank_test.ergm<-function(x,plot=FALSE){
  if(is.null(x$sample)) stop("Lack-of-fit test does not work for ERGMs that did not run MCMC.")

  eta <- ergm.eta(coef(x), x$etamap)
  MCMCeta <- ergm.eta(x$MCMCtheta, x$etamap)
  
  etasum <- c(x$sample %*% eta)
  w <- exp(c(x$sample %*% (eta-MCMCeta)))

  obs <- !is.null(x$sample.obs)
  if(obs){
    etasum.obs <- c(x$sample.obs %*% eta)
    w.obs <- exp(c(x$sample.obs %*% (eta-MCMCeta)))
  }
  
  if(plot){
    plot(density(etasum, weights=w),main=expression(paste("Density of ",log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y["obs"],";",hat(theta)))))),
         xlab=expression(log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y["obs"],";",hat(theta))))),zero.line=FALSE)
    #' @importFrom graphics abline
    if(!obs)
      abline(v=0)
    else
      plot(density(etasum.obs, weights=w.obs),zero.line=FALSE,add=TRUE)
  }

  if(!obs) sum(w*(etasum<0))/sum(w)
  else sum(mapply(function(eo,wo){
    wo*sum(w*(etasum<eo))/sum(w)
  }, etasum.obs, w.obs))/sum(w.obs)
}
