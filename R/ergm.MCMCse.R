ergm.MCMCse<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      model, 
                      lag.max=10, lag.max.miss=lag.max) {
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
  if(!is.null(statsmatrix.miss)){
   av.miss <- apply(statsmatrix.miss, 2, mean)
#  av.miss <- apply(statsmatrix.miss, 2, median)
   xsim.miss <- sweep(statsmatrix.miss, 2, av.miss,"-")
   xobs <- av.miss-av
  }
  theta.offset <- etamap$theta0
  theta.offset[!offsettheta] <- theta
#
# eta transformation
#
  eta0 <- ergm.eta(theta0, etamap)
  eta <-  ergm.eta(theta, etamap)
# etagrad <- ergm.etagrad(theta, etamap)
  etaparam <- eta-eta0
  etaparam <- etaparam[!offsetmap]
# etagrad <- etagrad[,!offsetmap,drop=FALSE]
# etagrad <- etagrad[!offsettheta,,drop=FALSE]
  xobs <- xobs[!offsetmap]
  xsim <- xsim[,!offsetmap, drop=FALSE]
  if(!is.null(statsmatrix.miss)){
   xsim.miss <- xsim.miss[,!offsetmap, drop=FALSE]
  }
#
# names(theta) <- dimnames(statsmatrix)[[2]]
  names(theta) <- names(theta0)
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
   E.miss <- 0
   lag.max.miss <- lag.max
   if(!is.null(statsmatrix.miss)){
#   z <- sweep(xsim.miss, 2, xobs, "-")
    z <- xsim.miss
    R <- acf(z, lag.max = lag.max.miss,
     type = "covariance", plot = FALSE)$acf
    if(dim(R)[2] > 1){
     part <- apply(R[-1,  ,  ,drop=FALSE], c(2, 3), sum)
    }else{
     part <- matrix(sum(R[-1,  ,  , drop=FALSE]))
    }
    cov.zbar.miss <- (R[1,  ,  ] + part + t(part))/nrow(xsim.miss)
    misspred <- xsim.miss %*% etaparam
    prob.miss <- max(misspred)
    prob.miss <- exp(misspred - prob.miss)
    prob.miss <- prob.miss/sum(prob.miss)
    E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
#   E.miss <- apply(xsim.miss,2,wtd.median,weight=prob.miss)
    htmp <- sweep(sweep(xsim.miss, 2, E.miss, "-"), 1, sqrt(prob.miss), "*")
#   htmp <- htmp %*% t(etagrad)
#   H.miss <- t(htmp) %*% htmp
    htmp.offset[,!offsetmap] <- htmp
    htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
    H.miss <- crossprod(htmp.offset, htmp.offset)
#   htmp <- crossprod(htmp, htmp)
#   H <- crossprod(t(etagrad),crossprod(htmp, t(etagrad)))
#   gradient.miss <- (xobs-E.miss) %*% t(etagrad)
    llg.offset[!offsetmap] <- xobs-E.miss
    gradient.miss <- ergm.etagradmult(theta.offset, llg.offset, etamap)
#   gradient.miss <- tcrossprod(xobs-E.miss, etagrad)
    cov.zbar.miss <- suppressWarnings(chol(cov.zbar.miss, pivot=TRUE))
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.miss
    cov.zbar.miss.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.miss.offset), etamap))
    cov.zbar.miss <- crossprod(cov.zbar.miss.offset, cov.zbar.miss.offset)
#   cov.zbar.miss <- as.matrix(etagrad %*% cov.zbar.miss %*% t(etagrad))
#   if(any(!offsettheta)){
#    H <- as.matrix(H[!offsettheta,])
#    H <- as.matrix(H[,!offsettheta])
#   }else{
#    H <- matrix(ncol=0,nrow=0)
#   }
#   cov.zbar <- cov.zbar[!offsettheta,,drop=FALSE]
#   cov.zbar <- cov.zbar[,!offsettheta,drop=FALSE]
    novar <- novar | (diag(H.miss)==0)
    H.miss <- H.miss[!novar,,drop=FALSE] 
    H.miss <- H.miss[,!novar,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[!novar,,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[,!novar,drop=FALSE] 
    gradient.miss <- gradient.miss[!novar] 
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
   if(!is.null(statsmatrix.miss)){
     gradient.full[!offsettheta][!novar] <- gradient - gradient.miss 
   }else{
     gradient.full[!offsettheta][!novar] <- gradient
   }
   mc.se <- rep(NA,length=length(theta))
   mc.se0 <- try(solve(H, cov.zbar), silent=TRUE)
   if(!(inherits(mc.se0,"try-error"))){
    mc.se0 <- try(diag(solve(H, t(mc.se0))), silent=TRUE)
    if(!(inherits(mc.se0,"try-error"))){
     if(!is.null(statsmatrix.miss)){
      mc.se.miss0 <- try(solve(H.miss, cov.zbar.miss), silent=TRUE)
      if(!(inherits(mc.se.miss0,"try-error"))){
       mc.se.miss0 <- try(diag(solve(H.miss, t(mc.se.miss0))), silent=TRUE)
       if(!inherits(mc.se.miss0,"try-error")){
        mc.se[!offsettheta][!novar] <- sqrt(mc.se0 + mc.se.miss0)
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
    if(!is.null(statsmatrix.miss)){
     test.hessian.miss <- try(any(is.na(sqrt(diag(robust.inverse(H.miss))))), silent=TRUE)
     if(inherits(test.hessian.miss,"try-error") || test.hessian.miss){
#                || detna(H.miss)< -25 ){
       hessian0 <- - H
     }else{
       hessian0 <-  H.miss-H
#      hessian0 <- -robust.inverse(var(xsim[,!novar,drop=FALSE]))-robust.inverse(var(xsim.miss[,!novar,drop=FALSE]))
#      hessian0 <- -var(xsim[,!novar])+var(xsim.miss[,!novar])
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
