#  File ergm/R/ergm.llik.R
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

llik.fun <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  # Convert theta to eta
  eta <- ergm.eta(theta.offset, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
  #
  basepred <- xsim %*% etaparam
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  llr <- sum(xobs * etaparam) - mb - varweight*vb
  #

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -800}

  # trustregion is the maximum value of llr that we actually trust.
  # So if llr>trustregion, return a value less than trustregion instead.
  if (llr>trustregion) {
    return(2*trustregion - llr)
  } else {
    return(llr)
  }
}



#llik.grad <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
#                      varweight=0.5, trustregion=20, eta0, etamap){
#  theta.offset <- etamap$init
#  theta.offset[!etamap$offsettheta] <- theta
#  eta <- ergm.eta(theta.offset, etamap)
#  etaparam <- eta-eta0
#  xsim[,etamap$offsetmap] <- 0
#  basepred <- xsim %*% etaparam
#  prob <- max(basepred)
#  prob <- probs*exp(basepred - prob)
#  prob <- prob/sum(prob)
#  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
#  llg <- xobs - E
#  llg[is.na(llg) | is.infinite(llg)] <- 0
##
## Penalize changes to trustregion
##
## llg <- llg - 2*(llg-trustregion)*(llg>trustregion)
## 
## The next lines are for the Hessian which optim does not use
##
## vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
## V <- -t(vtmp) %*% vtmp
## list(gradient=xobs-E,hessian=V)
## print(ergm.etagradmult(theta.offset, llg, etamap))
#  llg <- ergm.etagradmult(theta.offset, llg, etamap)
#  llg[!etamap$offsettheta]
#}



llik.grad <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                      varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
  basepred <- xsim %*% etaparam
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  llg <- xobs - E
  llg[is.na(llg) | is.infinite(llg)] <- 0
#
# Penalize changes to trustregion
#
# llg <- llg - 2*(llg-trustregion)*(llg>trustregion)
# 
# The next lines are for the Hessian which optim does not use
#
# vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# V <- -t(vtmp) %*% vtmp
# list(gradient=xobs-E,hessian=V)
# print(ergm.etagradmult(theta.offset, llg, etamap))
  llg.offset <- rep(0,length(etamap$offsetmap))
  llg.offset[!etamap$offsetmap] <- llg
  llg <- ergm.etagradmult(theta.offset, llg.offset, etamap)
# llg <- crossprod(llg, t(etagrad))
# llg <- etagrad %*% llg
# print(llg)
  llg[!etamap$offsettheta]
}




llik.hessian <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                         varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
# xsim[,etamap$offsettheta] <- 0
#
#    eta transformation
#
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
# etaparam <- etaparam[!etamap$offsettheta]
  basepred <- xsim %*% etaparam
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# htmp <- htmp %*% t(etagrad)
# H <- - t(htmp) %*% htmp
# htmp <- etagrad %*% t(htmp)
# H <- - htmp %*% t(htmp)
# htmp <- tcrossprod(etagrad, htmp)
  htmp.offset <- matrix(0, ncol = length(etamap$offsetmap), nrow = nrow(htmp))
# htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
# Notice the negative sign!
  H <- -crossprod(htmp.offset, htmp.offset)
# H <- -tcrossprod(htmp.offset, htmp.offset)
# H <- -tcrossprod(htmp, htmp)
# htmp <- tcrossprod(htmp, etagrad)
# H <- crossprod(htmp, htmp)
# H <- crossprod(t(etagrad),crossprod(H, t(etagrad)))
# He <- matrix(NA, ncol = length(etamap$offsettheta), 
#                  nrow = length(etamap$offsettheta))
# He[!etamap$offsettheta, !etamap$offsettheta] <- H
  He <- H[!etamap$offsettheta, !etamap$offsettheta, drop=FALSE]
  dimnames(He) <- list(names(theta), names(theta))
# H
  He
}





llik.hessian.naive <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                               varweight=0.5, eta0, etamap){
  xsim <- xsim[,!etamap$offsettheta, drop=FALSE]
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  etagrad <- ergm.etagrad(theta, etamap)

  # It does not matter if we subtract a constant from the (eta-eta0)^t g vector
  # prior to exponentiation.  So subtract to make maximum value = 0 for 
  # the sake of numberical stability:
  etaparam <- eta-eta0
  etaparam <- etaparam[!etamap$offsettheta]
  basepred <- xsim %*% etaparam
  w <- probs * exp(basepred - max(basepred))
  w <- w/sum(w)
  wtxsim <- sweep(xsim, 1, w, "*")
  swg <- colSums(wtxsim)
  H <- outer(swg,swg) - t(wtxsim) %*% xsim 
  
  # One last step, for the case of a curved EF:  Front- and back-multiply by
  # the gradient of eta(theta).
  H <- - crossprod(etagrad, crossprod(H, etagrad))
  He <- matrix(NA, ncol = length(theta), nrow = length(theta))
  He[!etamap$offsettheta, !etamap$offsettheta] <- H
  dimnames(He) <- list(names(theta), names(theta))
  He
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



llik.fun.EF <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
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
  llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}




llik.fun2 <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL, 
                      varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  llr <- sum(xobs * etaparam) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  llr
}





llik.grad2 <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                       varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  etaparam <- eta-eta0
  basepred <- xsim %*% etaparam
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- sum(xsim*as.vector(prob))
# 
#    The next lines are for the Hessian which optim does not use
# 
#    vtmp <- (xsim-E)*sqrt(prob)
#    V <- vtmp * vtmp
#    list(gradient=xobs-E,hessian=V)
  ergm.etagradmult(theta.offset, xobs-E, etamap)
}

llik.hessian2 <- llik.hessian





#####################################################################################
# --RETURNED--
#   llr: the naive log-likelihood ratio, l(eta) - l(eta0), approximated by
#        equation (5) of Hunter & Handcock (2010)
#####################################################################################

llik.fun3 <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL, 
                      varweight=0.5, trustregion=20, eta0, etamap){ # eqn (5)
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  deta <- matrix(eta-eta0,ncol=1) #px1
  basepred <- as.vector(xsim %*% deta) #nx1
  maxbase <- max(basepred)
  sum(xobs * deta) - maxbase - log(sum(probs*exp(basepred-maxbase)))
}



#####################################################################################
# --RETURNED--
#   the gradient of eta with ??  
#####################################################################################

llik.grad3 <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                       varweight=0.5, trustregion=20, eta0, etamap){ #eqn (11)
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  tmp <- probs * exp(basepred-maxbase)
  w <- tmp/sum(tmp)
  ergm.etagradmult (theta, xobs-apply(sweep(xsim,1,w,"*"),2,sum), etamap)
}



#####################################################################################
# --RETURNED--
#   He: the ?? Hessian matrix
#####################################################################################

llik.info3 <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                       varweight=0.5, eta0, etamap){ #eqn (12)
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta,etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  w <- probs * exp(basepred-maxbase) # not normalized; cov.wt will do that
  t(etagrad) %*% cov.wt(xsim,w,center=FALSE)$cov %*% etagrad
}


#####################################################################################
# --RETURNED--
#   the MCMC variance, as approximated by equation (13) of Hunter & Handcock (2010)
#   "sort of"
#####################################################################################

llik.mcmcvar3 <- function(theta, xobs, xsim, probs,  xsim.obs=NULL, probs.obs=NULL,
                          varweight=0.5, eta0, etamap){ #eqn (13) sort of
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  ebp <- exp(basepred-maxbase)
  sum(probs^2) * (sum(probs*ebp^2)/sum(probs*ebp)^2 - 1) 
# Assuming Z_i independent, eqn (13) becomes
# \sum_i p_i^2 \left( \frac{\sum_i p_i U_i^2}{[\sum_i p_i U_i]^2} - 1 \right)
}




#####################################################################################
# --RETURNED--
#   llr: the log-likelihood ratio of l(eta) - l(eta0) using ?? (what sort of approach)
#####################################################################################

llik.fun.median <- function(theta, xobs, xsim, probs, xsim.obs=NULL, probs.obs=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$init
  theta.offset[!etamap$offsettheta] <- theta
  # Convert theta to eta
  eta <- ergm.eta(theta.offset, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  # i.e., assuming that the network statistics are approximately normally 
  # distributed so that exp(eta * stats) is lognormal
  etaparam <- eta-eta0
# MSH: Is this robust?
  etaparam <- etaparam[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# The next line is right!
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
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
  llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * etaparam) - log(sum(probs*exp(xsim %*% etaparam)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}
