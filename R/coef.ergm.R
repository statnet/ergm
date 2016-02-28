#  File R/coef.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
coef.ergm <- function(object, ...){object$coef}
coefficients.ergm <- coef.ergm

vcov.ergm <- function(object, sources=c("all","model","estimation"), ...){
  sources <- match.arg(sources)

  src.mod <- sources %in% c("all", "model")
  src.est <- sources %in% c("all", "estimation")

  p <- length(object$coef)
  
  if(src.mod){
    if(is.null(object$hessian) && is.null(object$covar)){
      object$covar <- matrix(NA, p, p)
    }
    if(is.null(object$covar)){
      v.mod <- ginv(-object$hessian)
    }else{
      v.mod <- object$covar
    }
    v.mod[diag(v.mod)<0|is.infinite(object$coef),] <- NA
    v.mod[,diag(v.mod)<0|is.infinite(object$coef)] <- NA
    v.mod[object$offset,] <- 0
    v.mod[,object$offset] <- 0
     colnames(v.mod) <- rownames(v.mod) <- names(object$coef)
  }

  if(src.est){
    v.est<- if(is.null(object$est.cov)) matrix(0, p, p) else object$est.cov
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
