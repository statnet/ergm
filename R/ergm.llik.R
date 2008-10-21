#  File ergm/R/ergm.llik.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
llik.fun <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     penalty=0.5, trustregion=20, eta0, etamap){
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
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * x) - (mb + penalty*vb)
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
llik.grad <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                      penalty=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
  xsim[,etamap$offsetmap] <- 0
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
  llr <- xobs - E
  llr[is.na(llr) | is.infinite(llr)] <- 0
#
# Penalize changes to trustregion
#
  llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
# 
# The next lines are for the Hessian which optim does not use
#
# vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# V <- t(vtmp) %*% vtmp
# list(gradient=xobs-E,hessian=V)
  ergm.etagradmult(theta.offset, llr, etamap)
}

llik.hessian <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                         penalty=0.5, eta0, etamap){
# theta.offset <- etamap$theta0
# theta.offset[!etamap$offsettheta] <- theta
  namesx <- names(theta)
  xsim[,etamap$offsetmap] <- 0
#
#    eta transformation
#
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta, etamap)
  x <- eta-eta0
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
  htmp <- htmp %*% t(etagrad)
  H <- t(htmp) %*% htmp
  dimnames(H) <- list(namesx, namesx)
  -H
}

# DH:  This llik.exp function does not appear to be used anywhere.
# MSH: Yep, just another idea based on exp. family theory
# llik.exp <- function(theta, xobs, xsim, probs, 
#                      penalty=0.5, eta0, etamap){
#   eta <- ergm.eta(theta, etamap)
#   x <- eta-eta0
#   vb <- var(xsim)
#   llr <- -sum(xobs * x) + penalty*(t(x) %*% vb %*% x)
#   if(is.infinite(llr) | is.na(llr)){llr <- -800}
#   llr
# }

#
# Simple convergence
#
llik.fun2 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL, 
                      penalty=0.5, eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  x <- eta-eta0
  basepred <- xsim * x
  maxbase <- max(basepred)
  llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  llr
}

llik.grad2 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                       penalty=0.5, eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  x <- eta-eta0
  basepred <- xsim * x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- sum(xsim*prob)
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
                      penalty=0.5, eta0, etamap){ # eqn (5) 
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1) #px1
  basepred <- as.vector(xsim %*% deta) #nx1
  maxbase <- max(basepred)
  sum(xobs * deta) - maxbase - log(sum(probs*exp(basepred-maxbase)))
}

llik.grad3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                       penalty=0.5, eta0, etamap){ #eqn (11)
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  tmp <- probs * exp(basepred-maxbase)
  w <- tmp/sum(tmp)
  ergm.etagradmult (theta, xobs-apply(sweep(xsim,1,w,"*"),2,sum), etamap)
}


llik.info3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                       penalty=0.5, eta0, etamap){ #eqn (12)
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta,etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  w <- probs * exp(basepred-maxbase) # not normalized; cov.wt will do that
  t(etagrad) %*% cov.wt(xsim,w,center=FALSE)$cov %*% etagrad
}


llik.mcmcvar3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                          penalty=0.5, eta0, etamap){ #eqn (13) sort of
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1)
  basepred <- as.vector(xsim %*% deta)
  maxbase <- max(basepred)
  ebp <- exp(basepred-maxbase)
  sum(probs^2) * (sum(probs*ebp^2)/sum(probs*ebp)^2 - 1) 
# Assuming Z_i independent, eqn (13) becomes
# \sum_i p_i^2 \left( \frac{\sum_i p_i U_i^2}{[\sum_i p_i U_i]^2} - 1 \right)
}


