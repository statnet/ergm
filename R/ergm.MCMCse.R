#  File ergm/R/ergm.MCMCse.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#############################################################################
# The <ergm.MCMCse> function computes and returns the MCMC standard errors
################################################################################

ergm.MCMCse<-function(theta, init, statsmatrix, statsmatrix.obs,
                      model, 
                      lag.max=10, lag.max.obs=lag.max) {
  # Not sure why this is necessary, but:
  names(theta) <- names(init)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(init, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap
  etaparam <- (eta-eta0)[!offsetmap]

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

  # Take any theta offsets (values fixed at theta-1) into consideration
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

  cov.zbar.offset <- matrix(0, ncol = length(offsetmap), 
                            nrow = length(offsetmap))
  cov.zbar <- suppressWarnings(chol(cov.zbar, pivot=TRUE))
  cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar
  cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
  cov.zbar <- crossprod(cov.zbar.offset, cov.zbar.offset)

  # Identify canonical parameters corresponding to statistics that do not vary
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
    obspred <- xsim.obs %*% etaparam
    prob.obs <- exp(obspred - max(obspred))
    prob.obs <- prob.obs/sum(prob.obs)
    E.obs <- apply(sweep(xsim.obs, 1, prob.obs, "*"), 2, sum)
    htmp.obs <- sweep(sweep(xsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
    htmp.obs.offset <- matrix(0, ncol = length(offsetmap), nrow = nrow(htmp.obs))
    htmp.obs.offset[,!offsetmap] <- htmp.obs
    htmp.obs.offset <- t(ergm.etagradmult(theta.offset, t(htmp.obs.offset), etamap))
    H.obs <- crossprod(htmp.obs.offset, htmp.obs.offset)
    cov.zbar.obs.offset <- matrix(0, ncol = length(offsetmap), 
                                  nrow = length(offsetmap))
    cov.zbar.obs <- suppressWarnings(chol(cov.zbar.obs, pivot=TRUE))
    cov.zbar.obs.offset[!offsetmap,!offsetmap] <- cov.zbar.obs
    cov.zbar.obs.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.obs.offset), etamap))
    cov.zbar.obs <- crossprod(cov.zbar.obs.offset, cov.zbar.obs.offset)
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
  if(length(novar)==length(offsettheta)){
   novar <- novar | offsettheta
  }else{
   novar <- novar[!offsettheta]
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
  mc.cov <- matrix(NA,ncol=length(theta),nrow=length(theta))
  mc.cov0 <- try(solve(H, cov.zbar), silent=TRUE)
  if(!(inherits(mc.cov0,"try-error"))){
    mc.cov0 <- try(solve(H, t(mc.cov0)), silent=TRUE)
    if(!(inherits(mc.cov0,"try-error"))){
      if(!is.null(statsmatrix.obs)){
        mc.cov.obs0 <- try(solve(H.obs, cov.zbar.obs), silent=TRUE)
        if(!(inherits(mc.cov.obs0,"try-error"))){
          mc.cov.obs0 <- try(solve(H.obs, t(mc.cov.obs0)), silent=TRUE)
          if(!inherits(mc.cov.obs0,"try-error")){
            mc.cov[!novar,!novar] <- mc.cov0 + mc.cov.obs0
          }else{
            mc.cov[!novar,!novar] <- mc.cov0
          }
        }else{
          mc.cov[!novar,!novar] <- mc.cov0
        }
      }else{
        mc.cov[!novar,!novar] <- mc.cov0
      }
    }
  }
  colnames(mc.cov) <- names(theta)
  rownames(mc.cov) <- names(theta)
#
  return(list(mc.se=mc.se, mc.cov=mc.cov))
}
