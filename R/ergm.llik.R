#  File R/ergm.llik.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
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
#   xobs       : the vector of observed statistics 
#   xsim       : the matrix of simulated statistics 
#   probs      : the probability weight for each row of the stats matrix
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

llik.fun.lognormal <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, 
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Convert theta to eta
  eta <- ergm.eta(theta, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  llr <- sum(xobs * etaparam) - mb - varweight*vb
  #

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -800}

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
#   llg: the gradient of the log-likelihood using "naive" (importance sampling) method
#####################################################################################

llik.grad.IS <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, 
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + log(probs)
  
  # Calculate the estimating function values sans offset
  llg <- xobs - lweighted.mean(xsim, basepred)
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
llik.hessian.IS <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, 
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + log(probs)

  # Calculate the estimating function values sans offset
  esim <- t(ergm.etagradmult(theta, t(xsim), etamap))

  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H <- -lweighted.var(esim, basepred)

  dimnames(H) <- list(names(theta), names(theta))
  H
}


# DH:  This llik.exp function does not appear to be used anywhere.
# MSH: Yep, just another idea based on exp. family theory
# llik.exp <- function(theta, xobs, xsim, probs, 
#                      varweight=0.5, eta0, etamap){
#   eta <- ergm.eta(theta.offset, etamap)
#   etaparam <- eta-eta0
#   vb <- var(xsim)
#   llr <- -sum(xobs * etaparam) + varweight*(t(etaparam) %*% vb %*% etaparam)
#   if(is.infinite(llr) | is.na(llr)){llr <- -800}
#   llr
# }



#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ?? (what sort of approach)
#####################################################################################

llik.fun.EF <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
# The next line is right!
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# These lines standardize:
  basepred <- xsim %*% etaparam
#
  maxbase <- max(basepred)
  llr <- sum(xobs * etaparam) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (is.numeric(trustregion) && llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}




#####################################################################################
# --RETURNED--
#   llr: the "naive" log-likelihood ratio of l(eta) - l(eta0) using importance sampling (what sort of approach)
#            "Simple convergence"
#####################################################################################

llik.fun.IS <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL, 
                     varweight=0.5, trustregion=20,
                     dampening=FALSE,dampening.min.ess=100, dampening.level=0.1,
                     eta0, etamap){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights and the likelihood ratio
  basepred <- xsim %*% etaparam + log(probs)
  llr <- sum(xobs * etaparam) - log_sum_exp(basepred)
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
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ?? (what sort of approach)
#####################################################################################

llik.fun.median <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20,
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
# llr <- sum(xobs * etaparam) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=probs)
  sdb <- 1.4826*wtd.median(abs(basepred-mb), weight=probs)
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * etaparam) - (mb + varweight*sdb*sdb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  if (is.numeric(trustregion) && llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}

llik.fun.logtaylor <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL, 
	 	                     varweight=0.5, trustregion=20,  
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
	 	#   http://xianblog.wordpress.com/2010/09/24/effective-sample-size/ 
	 	    if(!is.na(ess) && {ess<dampening.min.ess}){ return(-Inf) } #.005*length(wts)) 
	 	  } 
	 	 
	 	  basepred <- xsim %*% etaparam 
	 	  ns <- length(basepred) 
	 	  mb <- sum(basepred*probs) 
	 	  vb <- sum(basepred*basepred*probs)-mb*mb 
	 	  skew <- sqrt(ns*(ns-1))*sum(((basepred-mb)^3)*probs)/(vb^(3/2)*(ns-2)) 
	 	  if(!is.finite(skew) | is.na(skew)){skew <- 0} 
	 	  part <- mb+vb/2 + sum(((basepred-mb)^3)*probs)/6 
	 	  llr <- sum(xobs * etaparam) - part 
	 	  # 
	 	 
	 	  # Simplistic error control;  -800 is effectively like -Inf: 
	 	  if(is.infinite(llr) | is.na(llr)){llr <- -800} 
	 	 
	 	  # trustregion is the maximum value of llr that we actually trust. 
	 	  # So if llr>trustregion, return a value less than trustregion instead. 
	 	  if (is.numeric(trustregion) && llr>trustregion) { 
	 	    return(2*trustregion - llr) 
	 	  } else { 
	 	    return(llr) 
	 	  } 
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
