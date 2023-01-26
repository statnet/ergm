#  File R/ergm.llik.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
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





#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using a lognormal
#        approximation; i.e., assuming that the network statistics are approximately
#        normally  distributed so that exp(eta * stats) is lognormal
#####################################################################################

llik.fun.lognormal <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Convert theta to eta
  eta <- ergm.eta(theta, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  mb <- lweighted.mean(basepred,lrowweights(xsim))
  vb <- lweighted.var(basepred,lrowweights(xsim))
  - mb - varweight*vb
}


#####################################################################################
# --RETURNED--
#   llg: the gradient of the log-likelihood using "naive" (importance sampling) method
#####################################################################################

llik.grad.IS <- function(theta, xsim,  xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)
  
  # Calculate the estimating function values sans offset
  llg <- - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta, llg, etamap))
  
  llg[is.na(llg)] <- 0 # Note: Before, infinite values would get zeroed as well. Let's see if this works.

  llg
}




#####################################################################################
# --RETURNED--
#   He: the naive approximation to the Hessian matrix - namely,
#          (sum_i w_i g_i)(sum_i w_i g_i)^t - sum_i(w_i g_i g_i^t),  where
#              g_i = the ith vector of statistics and
#              w_i = normalized version of exp((eta-eta0)^t g_i) so that sum_i w_i=1
#       this is equation (3.5) of Hunter & Handcock (2006)
#####################################################################################
llik.hessian.IS <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)

  # Calculate the estimating function values sans offset
  esim <- t(ergm.etagradmult(theta, t(xsim), etamap))

  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H <- -lweighted.var(esim, basepred)

  dimnames(H) <- list(names(theta), names(theta))
  H
}


#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ?? (what sort of approach)
#####################################################################################

llik.fun.EF <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  - maxbase - log(sum(rowweights(xsim)*exp(basepred-maxbase)))
}




#####################################################################################
# --RETURNED--
#   llr: the "naive" log-likelihood ratio of l(eta) - l(eta0) using importance sampling (what sort of approach)
#            "Simple convergence"
#####################################################################################

llik.fun.IS <- function(theta, xsim, xsim.obs=NULL, 
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights and the likelihood ratio
  basepred <- xsim %*% etaparam + lrowweights(xsim)
  - log_sum_exp(basepred)
}

#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ?? (what sort of approach)
#####################################################################################

llik.fun.median <- function(theta, xsim, xsim.obs=NULL,
                     varweight=0.5,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Convert theta to eta
  eta <- ergm.eta(theta, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  # i.e., assuming that the network statistics are approximately normally 
  # distributed so that exp(eta * stats) is lognormal
  etaparam <- eta-eta0
# These lines standardize:
  basepred <- xsim %*% etaparam
#
# maxbase <- max(basepred)
# llr <- - maxbase - log(sum(exp(lprobs)*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=rowweights(xsim))
  sdb <- 1.4826*wtd.median(abs(basepred-mb), weight=rowweights(xsim))
# 
# This is the log-likelihood ratio (and not its negative)
#
  - (mb + varweight*sdb*sdb)
}

llik.fun.logtaylor <- function(theta, xsim, xsim.obs=NULL, 
		                     varweight=0.5,
	 	                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1, 
	 	                     eta0, etamap){ 
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
