#  File ergm/R/ergm.llik.miss.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2010 the statnet development team
######################################################################
#
#   missing data code
#
llik.fun.miss <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.miss <- xsim.miss[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# The next line is right!
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
  misspred <- xsim.miss %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  mm <- sum(misspred*probs.miss)
  vm <- sum(misspred*misspred*probs.miss) - mm*mm
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
llik.grad.miss <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                      varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.miss <- xsim.miss[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  misspred <- xsim.miss %*% x
  prob.miss <- max(misspred)
  prob.miss <- probs.miss*exp(misspred - prob.miss)
  prob.miss <- prob.miss/sum(prob.miss)
  E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
  llg <- xobs + E.miss-E
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

llik.hessian.miss <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                         varweight=0.5, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  namesx <- names(theta)
#
#    eta transformation
#
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.miss <- xsim.miss[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  misspred <- xsim.miss %*% x
  prob.miss <- max(misspred)
  prob.miss <- probs.miss*exp(misspred - prob.miss)
  prob.miss <- prob.miss/sum(prob.miss)
  E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
  llg <- xobs + E.miss - E
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# htmp <- htmp %*% t(etagrad)
# H <- t(htmp) %*% htmp
  htmp.offset <- matrix(0, ncol = length(etamap$offsetmap), nrow = nrow(htmp))
  htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
  H <- crossprod(htmp.offset, htmp.offset)
##H <- etagrad %*% H %*% t(etagrad)
  htmp <- sweep(sweep(xsim.miss, 2, E.miss, "-"), 1, sqrt(prob.miss), "*")
# htmp <- htmp %*% t(etagrad)
# H.miss <- t(htmp) %*% htmp
##H.miss <- etagrad %*% H.miss %*% t(etagrad)
  htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
  H.miss <- crossprod(htmp.offset, htmp.offset)
  H <- H.miss-H
# He <- matrix(NA, ncol = length(theta), nrow = length(theta))
# He[!etamap$offsettheta, !etamap$offsettheta] <- H
# dimnames(He) <- list(namesx, namesx)
# He
  H
}
#
#  robust missing data code
#
llik.fun.miss.robust<- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.miss <- xsim.miss[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
  misspred <- xsim.miss %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=probs)
  vb <- 1.4826*wtd.median(abs(basepred-mb), weight=probs)
# print(c(mean(probs),mean(probs.miss),var(probs),var(probs.miss)))
  mm <- wtd.median(misspred, weight=probs.miss)
  vm <- 1.4826*wtd.median(abs(misspred-mm), weight=probs.miss)
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
llik.fun.miss.robust<- llik.fun.miss 
