#  File R/ergm.llik.obs.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#==========================================================================
# This file contains the following 14 functions for computing log likelihoods,
# gradients, hessians, and such for networks with observation process
#          <llik.fun.obs>         <llik.hessian.obs>
#          <llik.grad.obs>        <llik.fun.obs.robust>
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
#                offset coefficients; the not-offset terms are given by 'init'
#                of the 'etamap'
#   xobs       : the vector of observed statistics
#   xsim       : the matrix of simulated statistics
#   probs      : the probability weight for each row of the stats matrix
#   xsim.obs  : the 'xsim' counterpart for observation process
#   probs.obs : the 'probs' counterpart for observation process
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
llik.fun.obs.lognormal <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# The next line is right!
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# These lines standardize:
  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- sum(xobs * etaparam) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  mm <- sum(obspred*probs.obs)
  vm <- sum(obspred*obspred*probs.obs) - mm*mm
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * etaparam) + (mm + varweight*vm) - (mb + varweight*vb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (is.numeric(trustregion) && llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}



#####################################################################################
# --RETURNED--
#   llg: the gradient of the not-offset eta parameters with ??
#####################################################################################
llik.grad.obs.IS <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                      varweight=0.5, trustregion=20,
                      dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                      eta0, etamap){
  # Construct the parameter vector incl. offsets
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0
  
  # Calculate log-importance-weights (unconstrained)
  etaparam.no <- etaparam[!etamap$offsetmap]
  xsim.no <- xsim[,!etamap$offsetmap, drop=FALSE]
  basepred <- xsim.no %*% etaparam.no + log(probs)

  # Calculate log-importance-weights (constrained)
  xsim.obs.no <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  obspred <- xsim.obs.no %*% etaparam.no + log(probs.obs)
  
  llg <- xobs + lweighted.mean(xsim.obs, obspred) - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta.offset, llg, etamap))[,!etamap$offsettheta,drop=FALSE]

  llg[is.na(llg)] <- 0 # Note: Before, infinite values would get zeroed as well. Let's see if this works.
  
  llg
}



#####################################################################################
# --RETURNED--
#   He: the ?? Hessian matrix
#####################################################################################

llik.hessian.obs.IS <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Construct the parameter vector incl. offsets
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights (unconstrained)
  etaparam.no <- etaparam[!etamap$offsetmap]
  xsim.no <- xsim[,!etamap$offsetmap, drop=FALSE]
  basepred <- xsim.no %*% etaparam.no + log(probs)

  # Calculate log-importance-weights (constrained)
  xsim.obs.no <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  obspred <- xsim.obs.no %*% etaparam.no + log(probs.obs)

  # Calculate the estimating function values sans offset
  esim.no <- t(ergm.etagradmult(theta.offset, t(xsim), etamap))[,!etamap$offsettheta,drop=FALSE]
  osim.no <- t(ergm.etagradmult(theta.offset, t(xsim.obs), etamap))[,!etamap$offsettheta,drop=FALSE]
  
  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H.no <- lweighted.var(osim.no, obspred) - lweighted.var(esim.no, basepred)

  dimnames(H.no) <- list(names(theta), names(theta))
  H.no
}


#####################################################################################
# --RETURNED--
#   llr: the "naive" log-likelihood ratio of l(eta) - l(eta0) using importance sampling (what sort of approach)
#            "Simple convergence"
#####################################################################################

llik.fun.obs.IS <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL, 
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Construct the parameter vector incl. offsets
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0

  # Obtain canonical parameters incl. offsets and difference from sampled-from (unconstrained)
  etaparam.no <- etaparam[!etamap$offsetmap]
  xsim.no <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs.no <- xobs[!etamap$offsetmap]

    # Obtain canonical parameters incl. offsets and difference from sampled-from  (unconstrained)
  xsim.obs.no <- xsim[,!etamap$offsetmap, drop=FALSE]
  
  # Calculate log-importance-weights and the likelihood ratio
  basepred <- xsim %*% etaparam + log(probs)
  obspred <- xsim.obs %*% etaparam + log(probs.obs)
  llr <- sum(xobs.no * etaparam.no) + log_sum_exp(obspred) - log_sum_exp(basepred)
  
  # trustregion is the maximum value of llr that we actually trust.
  # So if llr>trustregion, return a value less than trustregion instead.
  if (is.numeric(trustregion) && llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}


#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ??
#                "robust obsing data code"
#####################################################################################

llik.fun.obs.robust<- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xsim.obs <- xsim.obs[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# These lines standardize:
  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- sum(xobs * etaparam) - maxbase - log(sum(probs*exp(basepred-maxbase)))
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
  llr <- sum(xobs * etaparam) + (mm + varweight*vm*vm) - (mb + varweight*vb*vb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (is.numeric(trustregion) && llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}

llik.fun.obs.robust<- llik.fun.obs.IS 
