#############################################################################
# The <ergm.MCMCse> function computes and returns the MCMC standard errors 
#
# --PARAMETERS--
#   theta           :  the vector of theta coefficients
#   theta0          :  the vector of initial theta coefficients
#   statsmatrix     :  the matrix of network statistics
#   statsmatrix.miss:  the matrix of network statistics on the network of
#                      missing edges
#   model           :  the model, as returned by <ergm.getmodel>
#   lag.max         :  the maximum lag at which to calculate the acf for the
#                      the network corresponding to 'statsmatrix'; default=10
#   lag.max.miss    :  the maximum lag at which to calculate the acf for the
#                      the network corresponding to 'statsmatrix.miss';
#                      default=lag.max
#
# --RETURNED--
#   list with one item:
#   mc.se: the vector of MCMC standard error estimates for each theta parameter
#
################################################################################

ergm.MCMCse<-function(theta, theta0, statsmatrix, statsmatrix.miss,
                      model, 
                      lag.max=10, lag.max.miss=lag.max) {
  # Not sure why this is necessary, but:
  names(theta) <- names(theta0)

  # Transform theta to eta
  etamap <- model$etamap
  eta0 <- ergm.eta(theta0, etamap)
  eta <-  ergm.eta(theta, etamap)
  offsettheta <- model$etamap$offsettheta
  offsetmap <- model$etamap$offsetmap
  etaparam <- (eta-eta0)[!offsetmap]

  # Center statmatrix (and statsmatrix.miss, if applicable)
  av <- apply(statsmatrix, 2, mean)
# av <- apply(statsmatrix,2,median)
  xsim <- sweep(statsmatrix, 2, av, "-")
  xobs <- -av
  if(!is.null(statsmatrix.miss)){
   av.miss <- apply(statsmatrix.miss, 2, mean)
#  av.miss <- apply(statsmatrix.miss, 2, median)
   xsim.miss <- sweep(statsmatrix.miss, 2, av.miss,"-")
   xsim.miss <- xsim.miss[,!offsetmap, drop=FALSE]
   xobs <- av.miss-av
  }
  xobs <- xobs[!offsetmap]
  xsim <- xsim[,!offsetmap, drop=FALSE]

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
    misspred <- xsim.miss %*% etaparam
    prob.miss <- exp(misspred - max(misspred))
    prob.miss <- prob.miss/sum(prob.miss)
    E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
    htmp <- sweep(sweep(xsim.miss, 2, E.miss, "-"), 1, sqrt(prob.miss), "*")
    htmp.offset[,!offsetmap] <- htmp
    htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
    H.miss <- crossprod(htmp.offset, htmp.offset)
    cov.zbar.miss <- suppressWarnings(chol(cov.zbar.miss, pivot=TRUE))
    cov.zbar.offset[!offsetmap,!offsetmap] <- cov.zbar.miss
    cov.zbar.offset <- t(ergm.etagradmult(theta.offset, t(cov.zbar.offset), etamap))
    cov.zbar.miss <- crossprod(cov.zbar.offset, cov.zbar.offset)
    novar <- novar | (diag(H.miss)==0)
    H.miss <- H.miss[!novar,,drop=FALSE] 
    H.miss <- H.miss[,!novar,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[!novar,,drop=FALSE] 
    cov.zbar.miss <- cov.zbar.miss[,!novar,drop=FALSE] 
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
#  return(list(mc.se=mc.se, mc.cov=NULL))
  return(list(mc.se=mc.se))
}
