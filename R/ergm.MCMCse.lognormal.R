ergm.MCMCse<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      H, H.miss, model, 
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
#
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
   cov.zbar.offset <- matrix(0, ncol = length(offsetmap), 
                                nrow = length(offsetmap))
   cov.zbar <- suppressWarnings(chol(cov.zbar, pivot=TRUE))
   cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar
   cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
   cov.zbar <- crossprod(cov.zbar.offset, cov.zbar.offset)
   novar <- diag(H)==0
#
#  Calculate the auto-covariance of the Conditional MCMC suff. stats.
#  and hence the Conditional MCMC s.e.
#
   E.miss <- 0
   lag.max.miss <- lag.max
   if(!is.null(statsmatrix.miss)){
    z <- xsim.miss
    R <- acf(z, lag.max = lag.max.miss,
     type = "covariance", plot = FALSE)$acf
    if(dim(R)[2] > 1){
     part <- apply(R[-1,  ,  ,drop=FALSE], c(2, 3), sum)
    }else{
     part <- matrix(sum(R[-1,  ,  , drop=FALSE]))
    }
    cov.zbar.miss <- (R[1,  ,  ] + part + t(part))/nrow(xsim.miss)
    cov.zbar.miss <- suppressWarnings(chol(cov.zbar.miss, pivot=TRUE))
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.miss
    cov.zbar.miss.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.miss.offset), etamap))
    cov.zbar.miss <- crossprod(cov.zbar.miss.offset, cov.zbar.miss.offset)
    novar <- novar | (diag(H.miss)==0)
    H.miss <- H.miss[!novar,,drop=FALSE] 
    H.miss <- H.miss[,!novar,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[!novar,,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[,!novar,drop=FALSE] 
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
    return(list(mc.se=mc.se, hessian=hessian))
   }
   cov.zbar <- cov.zbar[!novar,,drop=FALSE] 
   cov.zbar <- cov.zbar[,!novar,drop=FALSE] 
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
   list(mc.se=mc.se, hessian=hessian, covar=covar)
}
