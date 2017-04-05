#  File ergm/R/ergm.llik.obs.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###################################################################################
# Each of the <llik.X> functions computes either a likelihood function, a gradient
# function, or a Hessian matrix;  Each takes the same set of input parameters and
# these are described below; the return values differ however and so these are
# described above each function.   
####################################################################################
llik.fun.obs <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# The next line is right!
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
  obspred <- xsim.obs %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  mm <- sum(obspred*probs.obs)
  vm <- sum(obspred*obspred*probs.obs) - mm*mm
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * x) + (mm + varweight*vm) - (mb + varweight*vb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}




llik.grad.obs <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                      varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  obspred <- xsim.obs %*% x
  prob.obs <- max(obspred)
  prob.obs <- probs.obs*exp(obspred - prob.obs)
  prob.obs <- prob.obs/sum(prob.obs)
  E.obs <- apply(sweep(xsim.obs, 1, prob.obs, "*"), 2, sum)
  llg <- xobs + E.obs-E
  llg[is.na(llg) | is.infinite(llg)] <- 0
#
# Penalize changes to trustregion
#
# llg <- llg - 2*(llg-trustregion)*(llg>trustregion)
# 
# The next lines are for the Hessian which optim does not use
#
# vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# V <- t(vtmp) %*% vtmp
# list(gradient=xobs-E,hessian=V)
  llg.offset <- rep(0,length(etamap$offsetmap))
  llg.offset[!etamap$offsetmap] <- llg
  llg <- ergm.etagradmult(theta.offset, llg.offset, etamap)
  llg[!etamap$offsettheta]
}




llik.hessian.obs <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                         varweight=0.5, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
#
#    eta transformation
#
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)

  obspred <- xsim.obs %*% x
  prob.obs <- max(obspred)
  prob.obs <- probs.obs*exp(obspred - prob.obs)
  prob.obs <- prob.obs/sum(prob.obs)
  E.obs <- apply(sweep(xsim.obs, 1, prob.obs, "*"), 2, sum)

  llg <- xobs + E.obs - E
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# htmp <- htmp %*% t(etagrad)
# H <- t(htmp) %*% htmp
  htmp.offset <- matrix(0, ncol = length(etamap$offsetmap), nrow = nrow(htmp))
  htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
  H <- crossprod(htmp.offset, htmp.offset)
##H <- etagrad %*% H %*% t(etagrad)
  
  htmp.obs <- sweep(sweep(xsim.obs, 2, E.obs, "-"), 1, sqrt(prob.obs), "*")
# htmp <- htmp %*% t(etagrad)
# H.obs <- t(htmp) %*% htmp
##H.obs <- etagrad %*% H.obs %*% t(etagrad)
  htmp.obs.offset <- matrix(0, ncol = length(etamap$offsetmap), nrow = nrow(htmp.obs))
  htmp.obs.offset[,!etamap$offsetmap] <- htmp.obs
  htmp.obs.offset <- t(ergm.etagradmult(theta.offset, t(htmp.obs.offset), etamap))
  H.obs <- crossprod(htmp.obs.offset, htmp.obs.offset)
  
  He <- H[!etamap$offsettheta, !etamap$offsettheta, drop=FALSE]
  dimnames(He) <- list(names(theta), names(theta))
  He.obs <- H.obs[!etamap$offsettheta, !etamap$offsettheta, drop=FALSE]
  dimnames(He.obs) <- list(names(theta), names(theta))

  He.obs-He
# He <- matrix(NA, ncol = length(theta), nrow = length(theta))
# He[!etamap$offsettheta, !etamap$offsettheta] <- H
# dimnames(He) <- list(namestheta, namestheta)
# He
}





llik.fun.obs.robust<- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
  obspred <- xsim.obs %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=probs)
  vb <- 1.4826*wtd.median(abs(basepred-mb), weight=probs)
# print(c(mean(probs),mean(probs.obs),var(probs),var(probs.obs)))
  mm <- wtd.median(obspred, weight=probs.obs)
  vm <- 1.4826*wtd.median(abs(obspred-mm), weight=probs.obs)
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * x) + (mm + varweight*vm*vm) - (mb + varweight*vb*vb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}

llik.fun.obs.robust<- llik.fun.obs 
