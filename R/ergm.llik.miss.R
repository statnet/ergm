#==========================================================================
# This file contains the following 14 functions for computing log likelihoods,
# gradients, hessians, and such for networks with missing data
#          <llik.fun.miss>         <llik.hessian.miss>
#          <llik.grad.miss>        <llik.fun.miss.robust>
#=========================================================================





###################################################################################
# Each of the <llik.X> functions computes either a likelihood function, a gradient
# function, or a Hessian matrix;  Each takes the same set of input parameters and
# these are described below; the return values differ however and so these are
# described above each function.   
# 
#
# --PARAMETERS--
#   theta      : the vector of theta parameters; this is only used to solidify
#                offset coefficients; the not-offset terms are given by 'theta0'
#                of the 'etamap'
#   xobs       : the vector of observed statistics
#   xsim       : the matrix of simulated statistics
#   probs      : the probability weight for each row of the stats matrix
#   xsim.miss  : the 'xsim' counterpart for missing observations
#   probs.miss : the 'probs' counterpart for missing observations
#   varweight  : the weight by which the variance of the base predictions will be
#                scaled; the name of this param was changed from 'penalty' to better
#                reflect what this parameter actually is; default=0.5, which is the
#                "true"  weight, in the sense that the lognormal approximation is
#                given by
#                           sum(xobs * x) - mb - 0.5*vb
#   trustregion: the maximum value of the log-likelihood ratio that is trusted;
#                default=20
#   eta0       : the initial eta vector
#   etamap     : the theta -> eta mapping, as returned by <ergm.etamap>
#
#
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0)
#
####################################################################################




#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using a lognormal
#        approximation; i.e., assuming that the network statistics are approximately
#        normally  distributed so that exp(eta * stats) is lognormal
#####################################################################################                           
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



#####################################################################################
# --RETURNED--
#   llg: the gradient of the not-offset eta parameters with ??
#####################################################################################

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



#####################################################################################
# --RETURNED--
#   He: the ?? Hessian matrix
#####################################################################################

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





#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ??
#                "robust missing data code"
#####################################################################################

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
