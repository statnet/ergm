ergm.MCMCse<-function(theta, theta0, statsmatrix, statsmatrix.miss=NULL,
                      model, 
                      lag.max=50, lag.max.miss=lag.max) {
  av <- apply(statsmatrix, 2, mean)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
#
# eta transformation
#
  eta0 <- ergm.eta(theta0, model$etamap)
  eta <- ergm.eta(theta, model$etamap)
  etagrad <- ergm.etagrad(theta, model$etamap)
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
   detna <- function(x){x <- determinant(x); if(is.na(x[1])){x <- -40};x}
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
       mc.se[!novar] <- sqrt(mc.se0)    
   }
   names(mc.se) <- names(theta)
#
# use the exact Hessian if possible
# 
   test.hessian <- try(any(is.na(sqrt(diag(robust.inverse(V))))), silent=TRUE)
   if(inherits(test.hessian,"try-error") || test.hessian){
    hessian0 <- robust.inverse(var(xsim[,!novar]))
   }else{
     hessian0 <- - V    
   }
   hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
   hessian[!novar, !novar] <- hessian0
   dimnames(hessian) <- list(names(theta),names(theta))
   list(mc.se=mc.se, hessian=hessian, gradient=gradient)
}
