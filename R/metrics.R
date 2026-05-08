#  File R/ergm.llik.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#=================================================================================
# This file contains the following 14 functions for computing log likelihoods,
# gradients, hessians, and such:
#      <llik.fun>            <llik.fun.EF>     <llik.grad3>
#      <llik.grad>           <llik.fun2>       <llik.info3>
#      <llik.hessian>        <llik.grad2>      <llik.mcmcvar3>
#      <llik.hessian.naive>  <llik.hessian2>   <llik.fun.median>
#      <llik.exp>            <llik.fun3>
#=================================================================================




###################################################################################
# Each of the <llik.X> functions computes either a likelihood function, a gradient
# function, or a Hessian matrix;  Each takes the same set of input parameters and
# these are described below; the return values differ however and so these are
# described above each function.
#
# --PARAMETERS--
#   theta      : the vector of theta parameters; this is only used to solidify
#                offset coefficients; the not-offset terms are given by 'init'
#                of the 'etamap'
#   xsim       : the matrix of simulated statistics 
#   probs      : the probability weight for each row of the stats matrix
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
# --IGNORED PARAMETERS--
#   xsim.obs  : the 'xsim' counterpart for observation process; default=NULL
#   probs.obs : the 'probs' counterpart for observation process; default=NULL
#
####################################################################################

basepred <- function(theta, xsim, eta0, etamap, rowweights) {
  eta <- ergm.eta(theta, etamap)
  deta <- eta - eta0
  structure(xsim %*% deta + if (rowweights) lrowweights(xsim) else 0,
            deta = deta)
}


#####################################################################################
# Lognormal approximation:
# xsim ~ N ->
# dl ~~ -log(sum(exp(deta * xsim))) ~~ - mean(deta*xsim) - var(deta*xsim)/2
# l' ~~ - eta' * mean(xsim) - eta' * var(xsim) * deta
# l'' ~~ - eta' * var(xsim) * eta' = - var(xsim * eta')
# varweight is if we want to use something other than 1/2 for some reason.
#####################################################################################

llik.fun.lognormal <- function(theta, xsim, eta0, etamap, varweight = 0.5, ...) {
  basepred <- basepred(theta, xsim, eta0, etamap, FALSE)
  mb <- lweighted.mean(basepred,lrowweights(xsim))
  vb <- lweighted.var(basepred,lrowweights(xsim))
  - mb - varweight*vb
}

llik.grad.lognormal <- function(theta, xsim, eta0, etamap, varweight = 0.5, ...) {
  eta <- ergm.eta(theta, etamap)
  deta <- eta - eta0

  # TODO: These can be precomputed.
  m <- lweighted.mean(xsim, lrowweights(xsim))
  v <- lweighted.var(xsim, lrowweights(xsim))

  drop(m + 2 * varweight * v %*% deta) |>
    ergm.etagradmult(theta, v = _, etamap) |>
    t() |>
    replace(is.na, 0) |>
    (`-`)()
}

llik.hessian.lognormal <- function(theta, xsim, eta0, etamap, varweight = 0.5, ...) {
  ergm.etagradmultt(theta, xsim, etamap) |>
    lweighted.var(lrowweights(xsim)) |>
    (`*`)(-2 * varweight)
}


#####################################################################################
# "Naive" (Importance Sampling) approximation:
#####################################################################################

llik.grad.IS <- function(theta, xsim, eta0, etamap, ...) {
  basepred(theta, xsim, eta0, etamap, TRUE) |>
    lweighted.mean(xsim, logw = _) |>
    ergm.etagradmult(theta, v = _, etamap) |>
    t() |>
    replace(is.na, 0) |>
    (`-`)()
}

llik.hessian.IS <- function(theta, xsim, eta0, etamap, ...) {
  basepred <- basepred(theta, xsim, eta0, etamap, TRUE)
  ergm.etagradmultt(theta, xsim, etamap) |> lweighted.var(basepred) |> (`-`)()
}

llik.fun.IS <- function(theta, xsim, eta0, etamap, ...) {
  basepred(theta, xsim, eta0, etamap, TRUE) |>  log_sum_exp() |> (`-`)()
}

#####################################################################################
# Median: as Lognormal, but using median and MAD instead of mean and variance.
#####################################################################################

llik.fun.median <- function(theta, xsim, eta0, etamap, varweight = 0.5, ...) {
  basepred <- basepred(theta, xsim, eta0, etamap, FALSE)

  # alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=rowweights(xsim))
  sdb <- 1.4826*wtd.median(abs(basepred-mb), weight=rowweights(xsim))
# 
# This is the log-likelihood ratio (and not its negative)
#
  - (mb + varweight*sdb*sdb)
}

#####################################################################################
# LogTaylor: an experimental method, currently unclear what to do with
# it.
#####################################################################################

llik.fun.logtaylor <- function(theta, xsim, eta0, etamap,
                               dampening = FALSE, dampening.min.ess = 100, dampening.level = 0.1,
                               ...) {
	 	  # Convert theta to eta 
	 	  eta <- ergm.eta(theta, etamap) 
 	 
 	  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation 
	 	  etaparam <- eta-eta0 
	 	  # 
	 	  if (dampening) { 
	 	    #if theta_extended is (almost) outside convex hull don't trust theta 
	 	    eta_extended <- eta + etaparam*dampening.level 
	 	    expon_extended <- xsim %*% (eta_extended - eta0) 
	 	    wts <- exp(expon_extended) 
	 	    ess <- ceiling(sum(wts)^2/sum(wts^2)) 
	 	#   https://xianblog.wordpress.com/2010/09/24/effective-sample-size/ 
	 	    if(!is.na(ess) && {ess<dampening.min.ess}){ return(-Inf) } #.005*length(wts)) 
	 	  } 
	 	 
	 	  basepred <- xsim %*% etaparam 
	 	  ns <- length(basepred) 
	 	  mb <- sum(basepred*rowweights(xsim)) 
	 	  vb <- sum(basepred*basepred*rowweights(xsim))-mb*mb 
	 	  skew <- sqrt(ns*(ns-1))*sum(((basepred-mb)^3)*rowweights(xsim))/(vb^(3/2)*(ns-2)) 
	 	  if(!is.finite(skew) | is.na(skew)){skew <- 0} 
	 	  part <- mb+vb/2 + sum(((basepred-mb)^3)*rowweights(xsim))/6 
		  return(- part)
	 	} 

	 	ergm.llik.wins <- function(x,trim=.05, na.rm=TRUE) { 
	 	    if (trim == 0){return(x)} 
	 	    if ((trim < 0) | (trim>0.5) ) 
 	        stop("trimming must be reasonable") 
	 	    qtrim <- quantile(x,c(trim,.5, 1-trim),na.rm = na.rm) 
	 	    xbot <- qtrim[1] 
	 	    xtop <- qtrim[3] 
	 	    x[x < xbot] <- xbot 
  	    x[x > xtop] <- xtop 
	 	    return(x) 
	 	} 
