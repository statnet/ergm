ergm.MCMCse<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      model, 
                      lag.max=50, lag.max.miss=lag.max) {
# Adjust for offset
  statsmatrix <- statsmatrix[!model$etamap$offsettheta,!model$etamap$offsettheta]
#
  av <- apply(statsmatrix, 2, mean)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
  if(!is.null(statsmatrix.miss)){
   statsmatrix.miss <- statsmatrix.miss[!model$etamap$offsettheta,!model$etamap$offsettheta]
   av.miss <- apply(statsmatrix.miss, 2, mean)
   xsim.miss <- sweep(statsmatrix.miss, 2, av.miss,"-")
   dav <- av.miss-av
  }
#
# eta transformation
#
  eta0 <- ergm.eta(theta0, model$etamap)
  eta <- ergm.eta(theta, model$etamap)
  etagrad <- ergm.etagrad(theta, model$etamap)
  eta0 <- eta0[!model$etamap$offsettheta]
  eta <- eta[!model$etamap$offsettheta]
  etagrad <- etagrad[!model$etamap$offsettheta,!model$etamap$offsettheta]
  etaparam <- eta-eta0
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
   prob <- exp(xsim %*% etaparam)
   prob <- prob/sum(prob)
   E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
   vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
   V <- t(vtmp) %*% vtmp
   gradient <- (xobs-E) %*% t(etagrad)
   V <- etagrad %*% V %*% t(etagrad)
   cov.zbar <- etagrad %*% cov.zbar %*% t(etagrad)  
#
#  Calculate the auto-covariance of the Conditional MCMC suff. stats.
#  and hence the Conditional MCMC s.e.
#
   E.miss <- 0
   if(!is.null(statsmatrix.miss)){
    z <- sweep(xsim.miss, 2, xobs, "-")
    lag.max.miss <- min(round(sqrt(nrow(xsim.miss))),lag.max.miss)
    R <- acf(z, lag.max = lag.max.miss,
     type = "covariance", plot = FALSE)$acf
    if(dim(R)[2] > 1){
      part <- apply(R[-1,  ,  ], c(2, 3), sum)
    }else{
      part <- matrix(sum(R[-1,  ,  ]))
    }
    cov.zbar.miss <- (R[1,  ,  ] + part + t(part))/nrow(xsim.miss)
    prob.miss <- exp(xsim.miss %*% etaparam)
    prob.miss <- prob.miss/sum(prob.miss)
    E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
    vtmp <- sweep(sweep(xsim.miss, 2, E.miss, "-"), 1, sqrt(prob.miss), "*")
    V.miss <- t(vtmp) %*% vtmp
#
    gradient <- (dav+E.miss-E) %*% t(etagrad)
    V.miss <- etagrad %*% V.miss %*% t(etagrad)
    cov.zbar.miss <- etagrad %*% cov.zbar.miss %*% t(etagrad)  
   }
   detna <- function(x){x <- det(x); if(is.na(x)){x <- -40};x}
   novar <- diag(V)==0
   V <- V[!novar,,drop=FALSE] 
   V <- V[,!novar,drop=FALSE] 
   if(all(dim(V)==c(0,0))){
    hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
    mc.se <- rep(NA,length=length(theta))
    return(list(mc.se=mc.se, hessian=hessian, gradient=mc.se))
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
    hessian0 <- robust.inverse(var(xsim[,!novar]))
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
   hessian[!model$etamap$offsettheta,!model$etamap$offsettheta][!novar, !novar] <- hessian0
   dimnames(hessian) <- list(names(theta),names(theta))
   covar <- matrix(NA, ncol=length(theta), nrow=length(theta))
   covar[!model$etamap$offsettheta,!model$etamap$offsettheta][!novar, !novar] <- robust.inverse(-hessian0)
   dimnames(covar) <- list(names(theta),names(theta))
   list(mc.se=mc.se, hessian=hessian, gradient=gradient, covar=covar)
}
