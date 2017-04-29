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
llik.fun.obs.lognormal <- function(theta, xobs, xsim, lprobs, xsim.obs=NULL, lprobs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
# These lines standardize:
  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- sum(xobs * etaparam) - maxbase - log(sum(exp(lprobs)*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- sum(basepred*exp(lprobs))
  vb <- sum(basepred*basepred*exp(lprobs)) - mb*mb
  mm <- sum(obspred*exp(lprobs.obs))
  vm <- sum(obspred*obspred*exp(lprobs.obs)) - mm*mm
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
llik.grad.obs.IS <- function(theta, xobs, xsim, lprobs,  xsim.obs=NULL, lprobs.obs=NULL,
                      varweight=0.5, trustregion=20,
                      dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                      eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  
  # Calculate log-importance-weights (unconstrained)
  basepred <- xsim %*% etaparam + lprobs

  # Calculate log-importance-weights (constrained)
  obspred <- xsim.obs %*% etaparam + lprobs.obs
  
  llg <- xobs + lweighted.mean(xsim.obs, obspred) - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta, llg, etamap))

  llg[is.na(llg)] <- 0 # Note: Before, infinite values would get zeroed as well. Let's see if this works.
  
  llg
}



#####################################################################################
# --RETURNED--
#   He: the ?? Hessian matrix
#####################################################################################

llik.hessian.obs.IS <- function(theta, xobs, xsim, lprobs, xsim.obs=NULL, lprobs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights (unconstrained)
  basepred <- xsim %*% etaparam + lprobs

  # Calculate log-importance-weights (constrained)
  obspred <- xsim.obs %*% etaparam + lprobs.obs

  # Calculate the estimating function values sans offset
  esim <- t(ergm.etagradmult(theta, t(xsim), etamap))
  osim <- t(ergm.etagradmult(theta, t(xsim.obs), etamap))
  
  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H <- lweighted.var(osim, obspred) - lweighted.var(esim, basepred)

  dimnames(H) <- list(names(theta), names(theta))
  H
}


#####################################################################################
# --RETURNED--
#   llr: the "naive" log-likelihood ratio of l(eta) - l(eta0) using importance sampling (what sort of approach)
#            "Simple convergence"
#####################################################################################

llik.fun.obs.IS <- function(theta, xobs, xsim, lprobs, xsim.obs=NULL, lprobs.obs=NULL, 
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  
  # Calculate log-importance-weights and the likelihood ratio
  basepred <- xsim %*% etaparam + lprobs
  obspred <- xsim.obs %*% etaparam + lprobs.obs
  llr <- sum(xobs * etaparam) + log_sum_exp(obspred) - log_sum_exp(basepred)
  
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

llik.fun.obs.robust<- function(theta, xobs, xsim, lprobs, xsim.obs=NULL, lprobs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- sum(xobs * etaparam) - maxbase - log(sum(exp(lprobs)*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=exp(lprobs))
  vb <- 1.4826*wtd.median(abs(basepred-mb), weight=exp(lprobs))
# print(c(mean(exp(lprobs)),mean(exp(lprobs.obs)),var(exp(lprobs)),var(exp(lprobs.obs))))
  mm <- wtd.median(obspred, weight=exp(lprobs.obs))
  vm <- 1.4826*wtd.median(abs(obspred-mm), weight=exp(lprobs.obs))
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
