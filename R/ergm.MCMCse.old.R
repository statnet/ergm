#  File ergm/R/ergm.MCMCse.old.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
ergm.MCMCse.old<-function(theta, init, statsmatrix, statsmatrix.obs,
                      model, 
                      lag.max=10, lag.max.obs=lag.max) {
#
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap
  etamap <- model$etamap
#
# Adjust for any offset
  av <- apply(statsmatrix, 2, mean)
# av <- apply(statsmatrix,2,median)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
  if(!is.null(statsmatrix.obs)){
   av.obs <- apply(statsmatrix.obs, 2, mean)
#  av.obs <- apply(statsmatrix.obs, 2, median)
   xsim.obs <- sweep(statsmatrix.obs, 2, av.obs,"-")
   xobs <- av.obs-av
  }
  theta.offset <- etamap$init
  theta.offset[!offsettheta] <- theta
#
# eta transformation
#
  eta0 <- ergm.eta(init, etamap)
  eta <-  ergm.eta(theta, etamap)
# etagrad <- ergm.etagrad(theta, etamap)
  etaparam <- eta-eta0
  etaparam <- etaparam[!offsetmap]
# etagrad <- etagrad[,!offsetmap,drop=FALSE]
# etagrad <- etagrad[!offsettheta,,drop=FALSE]
  xobs <- xobs[!offsetmap]
  xsim <- xsim[,!offsetmap, drop=FALSE]
  if(!is.null(statsmatrix.obs)){
   xsim.obs <- xsim.obs[,!offsetmap, drop=FALSE]
  }
#
# names(theta) <- dimnames(statsmatrix)[[2]]
  names(theta) <- names(init)
#
#  Calculate the auto-covariance of the MCMC suff. stats.
#  and hence the MCMC s.e.
#
#  require("ts", quietly = TRUE, keep.source = FALSE)
   z <- sweep(xsim, 2, xobs, "-")
   lag.max <- min(round(sqrt(nrow(xsim))),lag.max)
   if(nrow(xsim) > 1000){
    lag.max <- round(15*(1+1000/nrow(xsim)))
   }
   R <- acf(z, lag.max = lag.max,
    type = "covariance", plot = FALSE)$acf
   if(dim(R)[2] > 1){
     part <- apply(R[-1,  ,  ,drop=FALSE], c(2, 3), sum)
   }else{
     part <- matrix(sum(R[-1,  ,  , drop=FALSE]))
   }
   cov.zbar <- (R[1,  ,  ] + part + t(part))/nrow(xsim)
   basepred <- xsim %*% etaparam
   prob <- max(basepred)
   prob <- exp(basepred - prob)
   prob <- prob/sum(prob)
   E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
#  E <- apply(xsim,2,wtd.median,weight=prob)
   htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
   htmp.offset <- matrix(0, ncol = length(offsetmap), nrow = nrow(htmp))
   htmp.offset[,!offsetmap] <- htmp
   htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
   H <- crossprod(htmp.offset, htmp.offset)
#  htmp <- htmp %*% t(etagrad)
#  H <- t(htmp) %*% htmp
#  htmp <- crossprod(htmp, htmp)
#  H <- crossprod(t(etagrad),crossprod(htmp, t(etagrad)))
#  gradient <- (xobs-E) %*% t(etagrad)
#  gradient <- tcrossprod(xobs-E, etagrad)
   llg.offset <- rep(0,length(offsetmap))
   llg.offset[!offsetmap] <- xobs-E
   gradient <- ergm.etagradmult(theta.offset, llg.offset, etamap)
   cov.zbar.offset <- matrix(0, ncol = length(offsetmap), 
                                nrow = length(offsetmap))
   cov.zbar <- suppressWarnings(chol(cov.zbar, pivot=TRUE))
   cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar
   cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
   cov.zbar <- crossprod(cov.zbar.offset, cov.zbar.offset)
#  cov.zbar <- as.matrix(etagrad %*% cov.zbar %*% t(etagrad))
#  if(any(!offsettheta)){
#   H <- as.matrix(H[!offsettheta,])
#   H <- as.matrix(H[,!offsettheta])
#  }else{
#   H <- matrix(ncol=0,nrow=0)
#  }
#  cov.zbar <- cov.zbar[!offsettheta,,drop=FALSE]
#  cov.zbar <- cov.zbar[,!offsettheta,drop=FALSE]
   novar <- diag(H)==0
#
#  Calculate the auto-covariance of the Conditional MCMC suff. stats.
#  and hence the Conditional MCMC s.e.
#
   E.obs <- 0
   lag.max.obs <- lag.max
   if(!is.null(statsmatrix.obs)){
#   z <- sweep(xsim.obs, 2, xobs, "-")
    z <- xsim.obs
    R <- acf(z, lag.max = lag.max.obs,
     type = "covariance", plot = FALSE)$acf
    if(dim(R)[2] > 1){
     part <- apply(R[-1,  ,  ,drop=FALSE], c(2, 3), sum)
    }else{
     part <- matrix(sum(R[-1,  ,  , drop=FALSE]))
    }
    cov.zbar.obs <- (R[1,  ,  ] + part + t(part))/nrow(xsim.obs)
    obspred <- xsim.obs %*% etaparam
    prob.obs <- max(obspred)
    prob.obs <- exp(obspred - prob.obs)
    prob.obs <- prob.obs/sum(prob.obs)
    E.obs <- apply(sweep(xsim.obs, 1, prob.obs, "*"), 2, sum)
#   E.obs <- apply(xsim.obs,2,wtd.median,weight=prob.obs)
    htmp <- sweep(sweep(xsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
#   htmp <- htmp %*% t(etagrad)
#   H.obs <- t(htmp) %*% htmp
    htmp.offset[,!offsetmap] <- htmp
    htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
    H.obs <- crossprod(htmp.offset, htmp.offset)
#   htmp <- crossprod(htmp, htmp)
#   H <- crossprod(t(etagrad),crossprod(htmp, t(etagrad)))
#   gradient.obs <- (xobs-E.obs) %*% t(etagrad)
    llg.offset[!offsetmap] <- xobs-E.obs
    gradient.obs <- ergm.etagradmult(theta.offset, llg.offset, etamap)
#   gradient.obs <- tcrossprod(xobs-E.obs, etagrad)
    cov.zbar.obs <- suppressWarnings(chol(cov.zbar.obs, pivot=TRUE))
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.obs
    cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
    cov.zbar.obs <- crossprod(cov.zbar.offset, cov.zbar.offset)
#   cov.zbar.obs <- as.matrix(etagrad %*% cov.zbar.obs %*% t(etagrad))
#   if(any(!offsettheta)){
#    H <- as.matrix(H[!offsettheta,])
#    H <- as.matrix(H[,!offsettheta])
#   }else{
#    H <- matrix(ncol=0,nrow=0)
#   }
#   cov.zbar <- cov.zbar[!offsettheta,,drop=FALSE]
#   cov.zbar <- cov.zbar[,!offsettheta,drop=FALSE]
    novar <- novar | (diag(H.obs)==0)
    H.obs <- H.obs[!novar,,drop=FALSE] 
    H.obs <- H.obs[,!novar,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[!novar,,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[,!novar,drop=FALSE] 
    gradient.obs <- gradient.obs[!novar] 
   }
   detna <- function(x){x <- det(x); if(is.na(x)){x <- -40};x}
   if(nrow(H)==1){
    H <- as.matrix(H[!novar,]) 
    H <- as.matrix(H[,!novar]) 
   }else{
    H <- H[!novar,,drop=FALSE] 
    H <- H[,!novar,drop=FALSE] 
   }
   if(all(dim(H)==c(0,0))){
    hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
    mc.se <- rep(NA,length=length(theta))
    return(list(mc.se=mc.se, hessian=hessian, gradient=mc.se))
   }
   cov.zbar <- cov.zbar[!novar,,drop=FALSE] 
   cov.zbar <- cov.zbar[,!novar,drop=FALSE] 
   gradient <- gradient[!novar] 
   gradient.full <- rep(NA,length=length(theta))
   if(!is.null(statsmatrix.obs)){
     gradient.full[!offsettheta][!novar] <- gradient - gradient.obs 
   }else{
     gradient.full[!offsettheta][!novar] <- gradient
   }
   mc.se <- rep(NA,length=length(theta))
   mc.se0 <- try(solve(H, cov.zbar), silent=TRUE)
   if(!(inherits(mc.se0,"try-error"))){
    mc.se0 <- try(diag(solve(H, t(mc.se0))), silent=TRUE)
    if(!(inherits(mc.se0,"try-error"))){
     if(!is.null(statsmatrix.obs)){
      mc.se.obs0 <- try(solve(H.obs, cov.zbar.obs), silent=TRUE)
      if(!(inherits(mc.se.obs0,"try-error"))){
       mc.se.obs0 <- try(diag(solve(H.obs, t(mc.se.obs0))), silent=TRUE)
       if(!inherits(mc.se.obs0,"try-error")){
        mc.se[!offsettheta][!novar] <- sqrt(mc.se0 + mc.se.obs0)
       }else{
        mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
       }
      }else{
       mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
      }
     }else{
      mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
     }
    }
   }
   names(mc.se) <- names(theta)
#
# use the exact Hessian if possible
# 
   test.hessian <- try(any(is.na(sqrt(diag(robust.inverse(H))))), silent=TRUE)
   if(inherits(test.hessian,"try-error") || test.hessian){
    hessian0 <- robust.inverse(var(xsim[,!novar,drop=FALSE]))
   }else{
    if(!is.null(statsmatrix.obs)){
     test.hessian.obs <- try(any(is.na(sqrt(diag(robust.inverse(H.obs))))), silent=TRUE)
     if(inherits(test.hessian.obs,"try-error") || test.hessian.obs){
#                || detna(H.obs)< -25 ){
       hessian0 <- - H
     }else{
       hessian0 <-  H.obs-H
#      hessian0 <- -robust.inverse(var(xsim[,!novar,drop=FALSE]))-robust.inverse(var(xsim.obs[,!novar,drop=FALSE]))
#      hessian0 <- -var(xsim[,!novar])+var(xsim.obs[,!novar])
     }
    }else{
     hessian0 <- - H
    }
   }
#  hessian0 <- - H
   hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
#  hessian <- hessian[!offsettheta,!offsettheta,drop=FALSE]
#  hessian[!novar, !novar] <- hessian0
#  hessian[!offsettheta,!offsettheta,drop=FALSE][!novar, !novar,drop=FALSE] <- hessian0
   if(nrow(H)==1){
    hessian[!novar, !novar] <- hessian0
   }else{
    hessian[!offsettheta,!offsettheta][!novar, !novar] <- hessian0
    hessian[!offsettheta,][novar,] <- NA
    hessian[,!offsettheta][,novar] <- NA
   }
   dimnames(hessian) <- list(names(theta),names(theta))
   covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
#  covar <- covar[!offsettheta,!offsettheta,drop=FALSE]
#  covar[!novar, !novar] <- robust.inverse(-hessian0)
   if(nrow(H)==1){
    covar[!novar, !novar] <- robust.inverse(-hessian0)
   }else{
    covar[!offsettheta,!offsettheta][!novar, !novar] <- robust.inverse(-hessian0)
    covar[!offsettheta,][novar,] <- NA
    covar[,!offsettheta][,novar] <- NA
   }
   dimnames(covar) <- list(names(theta),names(theta))
   list(mc.se=mc.se, hessian=hessian, gradient=gradient.full, covar=covar)
}
