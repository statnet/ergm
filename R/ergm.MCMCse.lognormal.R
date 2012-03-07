#  File ergm/R/ergm.MCMCse.lognormal.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#############################################################################
# The <ergm.MCMCse.lognormal> function computes and returns the MCMC lognormal
# standard errors 
################################################################################

ergm.MCMCse.lognormal<-function(theta, init, statsmatrix, statsmatrix.obs,
                      H, H.obs, model, 
                      lag.max=10, lag.max.obs=lag.max) {
  # Not sure why this is necessary, but:
  names(theta) <- names(init)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(init, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap

  # Center statmatrix (and statsmatrix.obs, if applicable)
  av <- apply(statsmatrix, 2, mean)
# av <- apply(statsmatrix,2,median)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
  if(!is.null(statsmatrix.obs)){
   av.obs <- apply(statsmatrix.obs, 2, mean)
#  av.obs <- apply(statsmatrix.obs, 2, median)
   xsim.obs <- sweep(statsmatrix.obs, 2, av.obs,"-")
   xsim.obs <- xsim.obs[,!offsetmap, drop=FALSE]
   xobs <- av.obs-av
  }
  xobs <- xobs[!offsetmap]
  xsim <- xsim[,!offsetmap, drop=FALSE]

  # Take any theta offsets (values fixed at init) into consideration
  theta.offset <- etamap$init
  theta.offset[!offsettheta] <- theta[!offsettheta]

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
  # Note that some care may be required here, as H and cov.zbar may not be
  # the same dimension in case of a curved EF model, in which case this 
  # is probably the wrong function to call!
  novar <- diag(H)==0

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  E.obs <- 0
  lag.max.obs <- lag.max
  if(!is.null(statsmatrix.obs)){
    z <- xsim.obs
    R <- acf(z, lag.max = lag.max.obs, type = "covariance", plot = FALSE)$acf
    if(dim(R)[2] > 1){
      part <- apply(R[-1,  ,  ,drop=FALSE], c(2, 3), sum)
    }else{
      part <- matrix(sum(R[-1,  ,  , drop=FALSE]))
    }
    cov.zbar.obs <- (R[1,  ,  ] + part + t(part))/nrow(xsim.obs)
    cov.zbar.obs <- suppressWarnings(chol(cov.zbar.obs, pivot=TRUE))
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.obs
    cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
    cov.zbar.obs <- crossprod(cov.zbar.offset, cov.zbar.offset)
    novar <- novar | (diag(H.obs)==0)
    H.obs <- H.obs[!novar,,drop=FALSE] 
    H.obs <- H.obs[,!novar,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[!novar,,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[,!novar,drop=FALSE] 
  }
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
    return(mc.se)
  }
  cov.zbar <- cov.zbar[!novar,,drop=FALSE] 
  cov.zbar <- cov.zbar[,!novar,drop=FALSE] 
  mc.se <- rep(NA,length=length(theta))
  mc.se0 <- try(solve(H, cov.zbar), silent=TRUE)
  if(length(novar)==length(offsettheta)){
   novar <- novar | offsettheta
  }else{
   novar <- novar[!offsettheta]
  }
  if(!(inherits(mc.se0,"try-error"))){
    mc.se0 <- try(diag(solve(H, t(mc.se0))), silent=TRUE)
    if(!(inherits(mc.se0,"try-error"))){
      if(!is.null(statsmatrix.obs)){
        mc.se.obs0 <- try(solve(H.obs, cov.zbar.obs), silent=TRUE)
        if(!(inherits(mc.se.obs0,"try-error"))){
          mc.se.obs0 <- try(diag(solve(H.obs, t(mc.se.obs0))), silent=TRUE)
          if(!inherits(mc.se.obs0,"try-error")){
            mc.se[!novar] <- sqrt(mc.se0 + mc.se.obs0)
          }else{
            mc.se[!novar] <- sqrt(mc.se0)
          }
        }else{
          mc.se[!novar] <- sqrt(mc.se0)
        }
      }else{
        mc.se[!novar] <- sqrt(mc.se0)
      }
    }
  }
  names(mc.se) <- names(theta)
  return(mc.se)
}
