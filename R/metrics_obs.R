#  File R/ergm.llik.obs.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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

baseobspred <- function(theta, xsim, xsim.obs, eta0, etamap, rowweights) {
  eta <- ergm.eta(theta, etamap)
  deta <- eta - eta0
  structure(list(base = xsim %*% deta + if (rowweights) lrowweights(xsim) else 0,
                 obs = xsim.obs %*% deta + if (rowweights) lrowweights(xsim.obs) else 0),
            deta = deta)
}


#####################################################################################
# Lognormal approximation
#####################################################################################                           
llik.fun.obs.lognormal <- function(theta, xsim, xsim.obs, eta0, etamap,
                                   varweight = 0.5, ...) {
  pred <- baseobspred(theta, xsim, xsim.obs, eta0, etamap, FALSE)

  mb <- lweighted.mean(pred$base, lrowweights(xsim))
  vb <- lweighted.var(pred$base, lrowweights(xsim))
  mo <- lweighted.mean(pred$obs, lrowweights(xsim.obs))
  vo <- lweighted.var(pred$obs, lrowweights(xsim.obs), 0)

  (mo + varweight*vo) - (mb + varweight*vb)
}

llik.grad.obs.lognormal <- function(theta, xsim, xsim.obs, eta0, etamap, varweight = 0.5, ...) {
  eta <- ergm.eta(theta, etamap)
  deta <- eta - eta0
  
  # TODO: These can be precomputed.
  mb <- lweighted.mean(xsim, lrowweights(xsim))
  vb <- lweighted.var(xsim, lrowweights(xsim))
  mo <- lweighted.mean(xsim.obs, lrowweights(xsim.obs))
  vo <- lweighted.var(xsim.obs, lrowweights(xsim.obs))

  drop((mo + 2 * varweight * vo %*% deta) -
       (mb + 2 * varweight * vb %*% deta)) |>
    ergm.etagradmult(theta, v = _, etamap) |>
    t() |>
    replace(is.na, 0)
}

llik.hessian.obs.lognormal <- function(theta, xsim, xsim.obs, eta0, etamap, varweight = 0.5, ...) {
  vb <- ergm.etagradmultt(theta, xsim, etamap) |>
    lweighted.var(lrowweights(xsim))
  vo <- ergm.etagradmultt(theta, xsim.obs, etamap) |>
    lweighted.var(lrowweights(xsim.obs))

    2 * varweight * (vo - vb)
}

#####################################################################################
# "Naive" (Importance Sampling) approximation
#####################################################################################
llik.grad.obs.IS <- function(theta, xsim, xsim.obs, eta0, etamap, ...) {
  pred <- baseobspred(theta, xsim, xsim.obs, eta0, etamap, TRUE)

  (lweighted.mean(xsim.obs, pred$obs) - lweighted.mean(xsim, pred$base)) |>
    ergm.etagradmult(theta, v = _, etamap) |>
    t() |>
    replace(is.na, 0)
}

llik.hessian.obs.IS <- function(theta, xsim, xsim.obs, eta0, etamap, ...) {
  pred <- baseobspred(theta, xsim, xsim.obs, eta0, etamap, TRUE)

  vb <- ergm.etagradmultt(theta, xsim, etamap) |>
    lweighted.var(pred$base)
  vo <- ergm.etagradmultt(theta, xsim.obs, etamap) |>
    lweighted.var(pred$obs, 0)
  
  vo - vb
}

llik.fun.obs.IS <- function(theta, xsim, xsim.obs, eta0, etamap, ...) {
  pred <- baseobspred(theta, xsim, xsim.obs, eta0, etamap, TRUE)
  log_sum_exp(pred$obs) - log_sum_exp(pred$base)
}


#####################################################################################
# "Median"
#####################################################################################

llik.fun.obs.median <- function(theta, xsim, xsim.obs, eta0, etamap,
                                varweight = 0.5, ...) {
  pred <- baseobspred(theta, xsim, xsim.obs, eta0, etamap, FALSE)

  # alternative based on log-normal approximation
  mb <- wtd.median(pred$base, weight=rowweights(xsim))
  vb <- 1.4826*wtd.median(abs(pred$base-mb), weight=rowweights(xsim))
  mo <- wtd.median(pred$obs, weight=rowweights(xsim.obs))
  vo <- 1.4826*wtd.median(abs(pred$obs-mo), weight=rowweights(xsim.obs))

  (mo + varweight*vo*vo) - (mb + varweight*vb*vb)
}
