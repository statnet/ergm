# Note:  Former "penalty" argument has been changed to "varweight" to better
# reflect what it actually is.  The default value of 0.5 is the "true" weight,
# in the sense that the lognormal approximation is given by
# sum(xobs * x) - mb - 0.5*vb
llik.fun <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  # Convert theta to eta
  eta <- ergm.eta(theta.offset, etamap)

  # Calculate approximation to l(eta) - l(eta0) using a lognormal approximation
  # i.e., assuming that the network statistics are approximately normally 
  # distributed so that exp(eta * stats) is lognormal
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
  #
  basepred <- xsim %*% x
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  llr <- sum(xobs * x) - mb - varweight*vb
  #
  maxbase <- max(basepred)
  llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))

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

#llik.grad <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
#                      varweight=0.5, trustregion=20, eta0, etamap){
#  theta.offset <- etamap$theta0
#  theta.offset[!etamap$offsettheta] <- theta
#  eta <- ergm.eta(theta.offset, etamap)
#  x <- eta-eta0
#  xsim[,etamap$offsetmap] <- 0
#  basepred <- xsim %*% x
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
llik.grad <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                      varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
  basepred <- xsim %*% x
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

llik.hessian <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                         varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  namesx <- names(theta)
# xsim[,etamap$offsettheta] <- 0
#
#    eta transformation
#
  eta <- ergm.eta(theta.offset, etamap)
# etagrad <- ergm.etagrad(theta.offset, etamap)
  x <- eta-eta0
# MSH: Is this robust?
  x <- x[!etamap$offsetmap]
  xsim <- xsim[,!etamap$offsetmap, drop=FALSE]
  xobs <- xobs[!etamap$offsetmap]
# etagrad <- etagrad[,!etamap$offsetmap,drop=FALSE]
# etagrad <- etagrad[!etamap$offsettheta,,drop=FALSE]
#
# x <- x[!etamap$offsettheta]
  basepred <- xsim %*% x
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
  htmp.offset[,!etamap$offsetmap] <- htmp
  htmp.offset <- t(ergm.etagradmult(theta.offset, t(htmp.offset), etamap))
# Notice the negative sign!
  H <- -crossprod(htmp.offset, htmp.offset)
# H <- -tcrossprod(htmp.offset, htmp.offset)
# H <- -tcrossprod(htmp, htmp)
# htmp <- tcrossprod(htmp, etagrad)
# H <- crossprod(htmp, htmp)
# H <- crossprod(t(etagrad),crossprod(H, t(etagrad)))
  He <- matrix(NA, ncol = length(etamap$offsettheta), 
                   nrow = length(etamap$offsettheta))
  He[!etamap$offsettheta, !etamap$offsettheta] <- H
  dimnames(He) <- list(names(namesx), names(namesx))
# H
  He
}

# Use the naive approximation to the Hessian matrix.  
# Namely, (sum_i w_i g_i)(sum_i w_i g_i)^t - sum_i(w_i g_i g_i^t)
# where g_i is the ith vector of statistics and
# w_i = normalized version of exp((eta-eta0)^t g_i) so that sum_i w_i=1
# This is equation (3.5) of Hunter and Handcock (2006)
# Currently, llik.hessian.naive is merely a rewrite of llik.hessian
llik.hessian.naive <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                               varweight=0.5, eta0, etamap){
  namesx <- names(theta)
  xsim <- xsim[,!etamap$offsettheta, drop=FALSE]
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta, etamap)

  # It does not matter if we subtract a constant from the (eta-eta0)^t g vector
  # prior to exponentiation.  So subtract to make maximum value = 0 for 
  # the sake of numberical stability:
  x <- eta-eta0
  x <- x[!etamap$offsettheta]
  basepred <- xsim %*% x
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
  dimnames(He) <- list(names(namesx), names(namesx))
  He
}



# DH:  This llik.exp function does not appear to be used anywhere.
# MSH: Yep, just another idea based on exp. family theory
# llik.exp <- function(theta, xobs, xsim, probs, 
#                      varweight=0.5, eta0, etamap){
#   eta <- ergm.eta(theta, etamap)
#   x <- eta-eta0
#   vb <- var(xsim)
#   llr <- -sum(xobs * x) + varweight*(t(x) %*% vb %*% x)
#   if(is.infinite(llr) | is.na(llr)){llr <- -800}
#   llr
# }
llik.fun.EF <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# The next line is right!
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
#
  maxbase <- max(basepred)
  llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}

#
# Simple convergence
#
llik.fun2 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL, 
                      varweight=0.5, trustregion=20, eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  x <- eta-eta0
  basepred <- xsim %*% x
  maxbase <- max(basepred)
  llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  llr
}

llik.grad2 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                       varweight=0.5, trustregion=20, eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  x <- eta-eta0
  basepred <- xsim %*% x
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
  ergm.etagradmult(theta, xobs-E, etamap)
}

llik.hessian2 <- llik.hessian

##### New stuff:  (Based on Hunter and Handcock)

llik.fun3 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL, 
                      varweight=0.5, trustregion=20, eta0, etamap){ # eqn (5) 
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1) #px1
  basepred <- as.vector(xsim %*% deta) #nx1
  maxbase <- max(basepred)
  sum(xobs * deta) - maxbase - log(sum(probs*exp(basepred-maxbase)))
}

llik.grad3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                       varweight=0.5, trustregion=20, eta0, etamap){ #eqn (11)
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  tmp <- probs * exp(basepred-maxbase)
  w <- tmp/sum(tmp)
  ergm.etagradmult (theta, xobs-apply(sweep(xsim,1,w,"*"),2,sum), etamap)
}


llik.info3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                       varweight=0.5, eta0, etamap){ #eqn (12)
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta,etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  w <- probs * exp(basepred-maxbase) # not normalized; cov.wt will do that
  t(etagrad) %*% cov.wt(xsim,w,center=FALSE)$cov %*% etagrad
}


llik.mcmcvar3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
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

llik.fun.median <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     varweight=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# The next line is right!
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- wtd.median(basepred, weight=probs)
  sdb <- 1.4826*wtd.median(abs(basepred-mb), weight=probs)
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * x) - (mb + varweight*sdb*sdb)
  if(is.infinite(llr) | is.na(llr)){llr <- -800}
#
# Penalize changes to trustregion
#
  llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
#
# cat(paste("max, log-lik",maxbase,llr,"\n"))
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# cat(paste("log-lik",llr,aaa,"\n"))
# aaa
  llr
}
