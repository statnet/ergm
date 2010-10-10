ergm.MCMCse.lognormal<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      H, H.miss, model, 
                      lag.max=10, lag.max.miss=lag.max) {
  # Not sure why this is necessary, but:
  names(theta) <- names(theta0)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(theta0, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap

  # Center statmatrix (and statsmatrix.miss, if applicable)
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

  # Take any theta offsets (values fixed at theta0) into consideration
  theta.offset <- etamap$theta0
  theta.offset[!offsettheta] <- theta

  #  Calculate the auto-covariance of the MCMC suff. stats.
  #  and hence the MCMC s.e.    
  z <- sweep(xsim, 2, xobs, "-")
  lag.max <- min(round(sqrt(nrow(xsim))),lag.max)
  if(nrow(xsim) > 1000){
    lag.max <- round(15*(1+1000/nrow(xsim)))
  }
  R <- acf(z, lag.max = lag.max, type = "covariance", plot = FALSE)$acf
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
  
  # Identify canonical parameters corresponding to statistics that do not vary
  # DRH:  Temporarily commenting out "novar" lines until they work for curved EF models
  #TEMP# novar <- diag(H)==0

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  E.miss <- 0
  lag.max.miss <- lag.max
  if(!is.null(statsmatrix.miss)){
    z <- xsim.miss
    R <- acf(z, lag.max = lag.max.miss, type = "covariance", plot = FALSE)$acf
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
    #TEMP# novar <- novar | (diag(H.miss)==0)
    #TEMP# H.miss <- H.miss[!novar,,drop=FALSE] 
    #TEMP# H.miss <- H.miss[,!novar,drop=FALSE] 
    #TEMP# cov.zbar.miss <- cov.zbar.miss[!novar,,drop=FALSE] 
    #TEMP# cov.zbar.miss <- cov.zbar.miss[,!novar,drop=FALSE] 
  }
  #TEMP# if(nrow(H)==1){
  #TEMP#   H <- as.matrix(H[!novar,]) 
  #TEMP#   H <- as.matrix(H[,!novar]) 
  #TEMP# }else{
  #TEMP#   H <- H[!novar,,drop=FALSE] 
  #TEMP#   H <- H[,!novar,drop=FALSE] 
  #TEMP# }
  if(all(dim(H)==c(0,0))){
    hessian <- matrix(NA, ncol=length(theta), nrow=length(theta))
    mc.se <- rep(NA,length=length(theta))
    return(list(mc.se=mc.se, hessian=hessian))
  }
  #TEMP# cov.zbar <- cov.zbar[!novar,,drop=FALSE] 
  #TEMP# cov.zbar <- cov.zbar[,!novar,drop=FALSE] 
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
            #TEMP# mc.se[!offsettheta][!novar] <- sqrt(mc.se0 + mc.se.miss0)
            mc.se[!offsettheta] <- sqrt(mc.se0 + mc.se.miss0)
          }else{
            #TEMP# mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
            mc.se[!offsettheta] <- sqrt(mc.se0)
          }
        }else{
          #TEMP# mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
          mc.se[!offsettheta] <- sqrt(mc.se0)
        }
      }else{
        #TEMP# mc.se[!offsettheta][!novar] <- sqrt(mc.se0)
        mc.se[!offsettheta] <- sqrt(mc.se0)
      }
    }
  }
  names(mc.se) <- names(theta)
  list(mc.se=mc.se)
}
