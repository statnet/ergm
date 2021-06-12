#  File R/ergm.llik.obs.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
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
#   xsim       : the matrix of simulated statistics
#   probs      : the probability weight for each row of the stats matrix
#   xsim.obs  : the 'xsim' counterpart for observation process
#   probs.obs : the 'probs' counterpart for observation process
#   varweight  : the weight by which the variance of the base predictions will be
#                scaled; the name of this param was changed from 'penalty' to better
#                reflect what this parameter actually is; default=0.5, which is the
#                "true"  weight, in the sense that the lognormal approximation is
#                given by
#                           - mb - 0.5*vb
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
llik.fun.obs.lognormal <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
# These lines standardize:
  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- - maxbase - log(sum(rowweights(xsim))*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- lweighted.mean(basepred,lrowweights(xsim))
  vb <- lweighted.var(basepred,lrowweights(xsim))
  mm <- lweighted.mean(obspred,lrowweights(xsim.obs))
  vm <- lweighted.var(obspred,lrowweights(xsim.obs))
# 
# This is the log-likelihood ratio (and not its negative)
#
  (mm + varweight*vm) - (mb + varweight*vb)
}



#####################################################################################
# --RETURNED--
#   llg: the gradient of the not-offset eta parameters with ??
#####################################################################################
llik.grad.obs.IS <- function(theta, xsim,  xsim.obs=NULL,
                      varweight=0.5,
                      dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                      eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  
  # Calculate log-importance-weights (unconstrained)
  basepred <- xsim %*% etaparam + lrowweights(xsim)

  # Calculate log-importance-weights (constrained)
  obspred <- xsim.obs %*% etaparam + lrowweights(xsim.obs)
  
  llg <- lweighted.mean(xsim.obs, obspred) - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta, llg, etamap))

  llg[is.na(llg)] <- 0 # Note: Before, infinite values would get zeroed as well. Let's see if this works.
  
  llg
}



#####################################################################################
# --RETURNED--
#   He: the ?? Hessian matrix
#####################################################################################

llik.hessian.obs.IS <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights (unconstrained)
  basepred <- xsim %*% etaparam + lrowweights(xsim)

  # Calculate log-importance-weights (constrained)
  obspred <- xsim.obs %*% etaparam + lrowweights(xsim.obs)

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

llik.fun.obs.IS <- function(theta, xsim, xsim.obs=NULL, 
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  
  # Calculate log-importance-weights and the likelihood ratio
  basepred <- xsim %*% etaparam + lrowweights(xsim)
  obspred <- xsim.obs %*% etaparam + lrowweights(xsim.obs)
  log_sum_exp(obspred) - log_sum_exp(basepred)
}


#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ??
#                "robust obsing data code"
#####################################################################################

llik.fun.obs.robust<- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
#
# maxbase <- max(basepred)
# llr <- - maxbase - log(sum(rowweights(xsim))*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=rowweights(xsim))
  vb <- 1.4826*wtd.median(abs(basepred-mb), weight=rowweights(xsim))
# print(c(mean(rowweights(xsim))),mean(rowweights(xsim.obs)),var(rowweights(xsim))),var(rowweights(xsim.obs)))
  mm <- wtd.median(obspred, weight=rowweights(xsim.obs))
  vm <- 1.4826*wtd.median(abs(obspred-mm), weight=rowweights(xsim.obs))
# 
# This is the log-likelihood ratio (and not its negative)
#
  (mm + varweight*vm*vm) - (mb + varweight*vb*vb)
}

llik.fun.obs.robust<- llik.fun.obs.IS 
