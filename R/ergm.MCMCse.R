#  File R/ergm.MCMCse.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#############################################################################
# The <ergm.MCMCse> function computes and returns the MCMC standard errors
#
# --PARAMETERS--
#   theta           :  the vector of theta coefficients
#   init          :  the vector of initial theta coefficients
#   statsmatrix     :  the matrix of network statistics
#   statsmatrix.obs :  the matrix of network statistics on the constrained network
#   model           :  the model, as returned by <ergm.getmodel>
#
# --RETURNED--
#   the variance of the MCMC sampling as a list containing:
#     mc.se : the vector of MCMC standard error estimates for each theta parameter
#     mc.cov: the MCMC covariance matrix of the theta parameters
#
################################################################################

ergm.MCMCse<-function(theta, init, statsmatrix, statsmatrix.obs,
                      model) {
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

  # Take any theta offsets (values fixed at init) into consideration
  theta.offset <- etamap$init
  theta.offset[!offsettheta] <- theta[!offsettheta]

  # Calculate prediction probabilities that had been used in the last MCMLE update.
  basepred <- xsim %*% etaparam
  prob <- max(basepred)
  prob <- exp(basepred - prob)
  prob <- prob/sum(prob)

  # Estimate the Hessian.
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  #  E <- apply(xsim,2,wtd.median,weight=prob)
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
  htmp.offset <- matrix(0, ncol = length(offsetmap), nrow = nrow(htmp))
  htmp.offset[,!offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
  H <- crossprod(htmp.offset, htmp.offset)

  #  Calculate the auto-covariance of the MCMC suff. stats.
  #  and hence the MCMC s.e.
  z <- sweep(xsim, 2, xobs, "-")
  cov.zbar <- .ergm.mvar.spec0(z) * sum(prob^2)
  imp.factor <- mean(prob^2)
  cov.zbar.offset <- matrix(0, ncol = length(offsetmap), 
                            nrow = length(offsetmap))
  cov.zbar <- suppressWarnings(chol(cov.zbar, pivot=TRUE))
  pivot <- order(attr(cov.zbar, "pivot"))
  cov.zbar <-cov.zbar[, pivot]
  cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar
  cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
  cov.zbar <- crossprod(cov.zbar.offset, cov.zbar.offset)
  
  # Identify canonical parameters corresponding to statistics that do not vary
  novar <- diag(H)==0
  novar.offset <- rep(TRUE, length(offsettheta))
  novar.offset[!offsettheta] <- novar # Note that novar.offset == TRUE where offsettheta==TRUE as well.

  #  Calculate the auto-covariance of the Conditional MCMC suff. stats.
  #  and hence the Conditional MCMC s.e.
  E.obs <- 0
  if(!is.null(statsmatrix.obs)){
    obspred <- xsim.obs %*% etaparam
    prob.obs <- exp(obspred - max(obspred))
    prob.obs <- prob.obs/sum(prob.obs)
    E.obs <- apply(sweep(xsim.obs, 1, prob.obs, "*"), 2, sum)
    htmp.obs <- sweep(sweep(xsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
    htmp.obs.offset <- matrix(0, ncol = length(offsetmap), nrow = nrow(htmp.obs))
    htmp.obs.offset[,!offsetmap] <- htmp.obs
    htmp.obs.offset <- t(ergm.etagradmult(theta.offset, t(htmp.obs.offset), etamap))
    H.obs <- crossprod(htmp.obs.offset, htmp.obs.offset)

    z <- xsim.obs
    cov.zbar.obs <- .ergm.mvar.spec0(z) * sum(prob.obs^2)
    imp.factor <- mean(imp.factor, mean(prob.obs^2))
    cov.zbar.obs.offset <- matrix(0, ncol = length(offsetmap), 
                                  nrow = length(offsetmap))
    cov.zbar.obs <- suppressWarnings(chol(cov.zbar.obs, pivot=TRUE))
    pivot <- order(attr(cov.zbar.obs, "pivot"))
    cov.zbar.obs <-cov.zbar.obs[, pivot]
    cov.zbar.obs.offset[!offsetmap,!offsetmap] <- cov.zbar.obs
    cov.zbar.obs.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.obs.offset), etamap))
    cov.zbar.obs <- crossprod(cov.zbar.obs.offset, cov.zbar.obs.offset)
    novar <- novar | (diag(H.obs)==0)
    H.obs <- H.obs[!novar,,drop=FALSE] 
    H.obs <- H.obs[,!novar,drop=FALSE] 
    cov.zbar.obs <- cov.zbar.obs[!(novar|offsettheta),!(novar|offsettheta),drop=FALSE]
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
  cov.zbar <- cov.zbar[!(novar|offsettheta),!(novar|offsettheta),drop=FALSE]
  if(length(novar)==length(offsettheta)){
   novar <- novar | offsettheta
  }else{
   novar <- novar[!offsettheta]
  }

  if(inherits(try(solve(H)),"try-error")) warning("Approximate Hessian matrix is singular. Standard errors due to MCMC approximation of the likelihood cannot be evaluated. This is likely due to highly correlated model terms.")

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

  attr(mc.cov, "imp.factor") <- imp.factor
  
  mc.cov
}
