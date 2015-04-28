#  File R/ergm.MCMCse.lognormal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
#############################################################################
# The <ergm.MCMCse.lognormal> function computes and returns the MCMC lognormal
# standard errors 
#
# --PARAMETERS--
#   theta           :  the vector of theta coefficients
#   init          :  the vector of initial theta coefficients
#   statsmatrix     :  the matrix of network statistics
#   statsmatrix.obs :  the matrix of network statistics on the constrained network
#   H               :  the Hessian matrix
#   H.obs           :  the Hessian matrix on the constrained network
#   model           :  the model, as returned by <ergm.getmodel>
#
# --RETURNED--
#   mc.se: the vector of MCMC lognormal standard error estimates for each theta
#          parameter
#
################################################################################

ergm.MCMCse.lognormal<-function(theta, init, statsmatrix, statsmatrix.obs,
                      H, H.obs, model) {
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
  cov.zbar <- .ergm.mvar.spec0(z) / nrow(z)
  cov.zbar.offset <- matrix(0, ncol = length(offsetmap), 
                            nrow = length(offsetmap))
  cov.zbar <- suppressWarnings(chol(cov.zbar,pivot=TRUE))
  pivot <- order(attr(cov.zbar, "pivot"))
  cov.zbar <-cov.zbar[, pivot]
  cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar
  cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
  cov.zbar <- crossprod(cov.zbar.offset, cov.zbar.offset)

  # Identify canonical parameters corresponding to statistics that do not vary
  # Note that some care may be required here, as H and cov.zbar may not be
  # the same dimension in case of a curved EF model, in which case this 
  # is probably the wrong function to call!
  novar <- diag(H)==0
  novar.offset <- rep(TRUE, length(offsettheta))
  novar.offset[!offsettheta] <- novar # Note that novar.offset == TRUE where offsettheta==TRUE as well.

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  E.obs <- 0
  if(!is.null(statsmatrix.obs)){
    z <- xsim.obs
    cov.zbar.obs <- .ergm.mvar.spec0(z) / nrow(z)
    cov.zbar.obs <- suppressWarnings(chol(cov.zbar.obs, pivot=TRUE))
    pivot <- order(attr(cov.zbar.obs, "pivot"))
    cov.zbar.obs <-cov.zbar.obs[, pivot]
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.obs
    cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
    cov.zbar.obs <- crossprod(cov.zbar.offset, cov.zbar.offset)

    novar.obs <- diag(H.obs)==0
    novar.offset.obs <- rep(TRUE, length(offsettheta))
    novar.offset.obs[!offsettheta] <- novar.obs

    novar.offset <- novar.offset | novar.offset.obs
    novar <- novar | novar.obs
    
    H.obs <- H.obs[!novar,,drop=FALSE] 
    H.obs <- H.obs[,!novar,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[!(novar.offset),!(novar.offset),drop=FALSE] 
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
    return(matrix(NA, length(theta), length(theta)))
  }
  cov.zbar <- cov.zbar[!(novar.offset),!(novar.offset),drop=FALSE]

  mc.cov <- matrix(NA,ncol=length(theta),nrow=length(theta))

  if(is.null(statsmatrix.obs)){
    mc.cov0 <- try(solve(H, cov.zbar), silent=TRUE)
    if(!(inherits(mc.cov0,"try-error"))){
      mc.cov0 <- try(solve(H, t(mc.cov0)), silent=TRUE)
      if(!(inherits(mc.cov0,"try-error"))){
        mc.cov[!novar.offset,!novar.offset] <- mc.cov0
      }
    }
  }else{
    H <- H.obs - H # Bread^-1
    cov.zbar <- cov.zbar + cov.zbar.obs # Filling
    
    mc.cov0 <- try(solve(H, cov.zbar), silent=TRUE)
    if(!(inherits(mc.cov0,"try-error"))){
      mc.cov0 <- try(solve(H, t(mc.cov0)), silent=TRUE)
      if(!(inherits(mc.cov0,"try-error"))){
        mc.cov[!novar.offset,!novar.offset] <- mc.cov0
      }
    }
  }
  colnames(mc.cov) <- names(theta)
  rownames(mc.cov) <- names(theta)

  mc.cov
}
