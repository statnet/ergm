ergm.MCMCse.lognormal<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      V, V.miss, model, 
                      lag.max=10, lag.max.miss=lag.max) {
# Adjust for any offset
  av <- apply(statsmatrix,2,median)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
  if(!is.null(statsmatrix.miss)){
   av.miss <- apply(statsmatrix.miss, 2, median)
   xsim.miss <- sweep(statsmatrix.miss, 2, av.miss,"-")
   xobs <- av.miss-av
  }
#
# eta transformation
#
  eta0 <- ergm.eta(theta0, model$etamap)
  eta <- ergm.eta(theta, model$etamap)
  etagrad <- ergm.etagrad(theta, model$etamap)
  etaparam <- eta-eta0
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
   cov.zbar <- as.matrix(etagrad %*% cov.zbar %*% t(etagrad))
   cov.zbar <- cov.zbar[!model$etamap$offsettheta,,drop=FALSE]
   cov.zbar <- cov.zbar[,!model$etamap$offsettheta,drop=FALSE]
   novar <- diag(V)==0
#
#  Calculate the auto-covariance of the Conditional MCMC suff. stats.
#  and hence the Conditional MCMC s.e.
#
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
    cov.zbar.miss <- etagrad %*% cov.zbar.miss %*% t(etagrad)  
    cov.zbar.miss <- cov.zbar.miss[,!model$etamap$offsettheta,drop=FALSE]
    cov.zbar.miss <- cov.zbar.miss[!model$etamap$offsettheta,,drop=FALSE]
    novar <- novar & diag(V.miss)==0
    V.miss <- V.miss[!novar,,drop=FALSE] 
    V.miss <- V.miss[,!novar,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[!novar,,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[,!novar,drop=FALSE] 
   }
   detna <- function(x){x <- det(x); if(is.na(x)){x <- -40};x}
   if(nrow(V)==1){
    V <- as.matrix(V[!novar,]) 
    V <- as.matrix(V[,!novar]) 
   }else{
    V <- V[!novar,,drop=FALSE] 
    V <- V[,!novar,drop=FALSE] 
   }
   if(all(dim(V)==c(0,0))){
    hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
    mc.se <- rep(NA,length=length(theta))
    return(list(mc.se=mc.se, hessian=hessian))
   }
   cov.zbar <- cov.zbar[!novar,,drop=FALSE] 
   cov.zbar <- cov.zbar[,!novar,drop=FALSE] 
   mc.se <- rep(NA,length=length(theta))
   mc.se0 <- try(diag(solve(V, t(solve(V, cov.zbar)))), silent=TRUE)
   if(!(inherits(mc.se0,"try-error") || detna(V)< -20)){
    if(!is.null(statsmatrix.miss)){
      mc.se.miss0 <- try(diag(solve(V.miss, t(solve(V.miss, cov.zbar.miss)))),
                         silent=TRUE)
      if(inherits(mc.se.miss0,"try-error") || detna(V.miss)< -20){
       mc.se[!model$etamap$offsettheta][!novar] <- sqrt(mc.se0)
      }else{
       mc.se[!model$etamap$offsettheta][!novar] <- sqrt(mc.se0 + mc.se.miss0)
      }
    }else{
       mc.se[!model$etamap$offsettheta][!novar] <- sqrt(mc.se0)
    }
   }
   names(mc.se) <- names(theta)
#
# use the exact Hessian if possible
# 
   test.hessian <- try(any(is.na(sqrt(diag(robust.inverse(V))))), silent=TRUE)
   if(inherits(test.hessian,"try-error") || test.hessian){
    hessian0 <- robust.inverse(var(xsim[,!novar,drop=FALSE]))
   }else{
    if(!is.null(statsmatrix.miss)){
     test.hessian.miss <- try(any(is.na(sqrt(diag(robust.inverse(V.miss))))), silent=TRUE)
     if(inherits(test.hessian.miss,"try-error") || test.hessian.miss
                 || detna(V.miss)< -20 ){
       hessian0 <- - V
     }else{
       hessian0 <-  V.miss-V
     }
    }else{
     hessian0 <- - V
    }
   }
   hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
   if(nrow(V)==1){
    hessian[!novar, !novar] <- hessian0
   }else{
    hessian[!model$etamap$offsettheta,!model$etamap$offsettheta][!novar, !novar] <- hessian0
    hessian[!model$etamap$offsettheta,][novar,] <- NA
    hessian[,!model$etamap$offsettheta][,novar] <- NA
   }
   dimnames(hessian) <- list(names(theta),names(theta))
   covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
   if(nrow(V)==1){
    covar[!novar, !novar] <- robust.inverse(-hessian0)
   }else{
    covar[!model$etamap$offsettheta,!model$etamap$offsettheta][!novar, !novar] <- robust.inverse(-hessian0)
    covar[!model$etamap$offsettheta,][novar,] <- NA
    covar[,!model$etamap$offsettheta][,novar] <- NA
   }
   dimnames(covar) <- list(names(theta),names(theta))
   list(mc.se=mc.se, hessian=hessian, covar=covar)
}
