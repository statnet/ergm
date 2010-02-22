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
  basepred <- xsim %*% x
  mb <- sum(basepred*probs)
  vb <- sum(basepred*basepred*probs) - mb*mb
  llr <- sum(xobs * x) - mb - varweight*vb

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

llik.grad <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                      varweight=0.5, trustregion=20, eta0, etamap){
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
# llr <- llr - 2*(llr-trustregion)*(llr>trustregion)
# 
# The next lines are for the Hessian which optim does not use
#
# vtmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
# V <- t(vtmp) %*% vtmp
# list(gradient=xobs-E,hessian=V)
# print(ergm.etagradmult(theta.offset, llr, etamap))
  llr <- ergm.etagradmult(theta.offset, llr, etamap)
  llr[!etamap$offsetmap]
}

llik.hessian <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                         varweight=0.5, eta0, etamap){
# theta.offset <- etamap$theta0
# theta.offset[!etamap$offsettheta] <- theta
  namesx <- names(theta)
# xsim[,etamap$offsettheta] <- 0
  xsim <- xsim[,!etamap$offsettheta]
#
#    eta transformation
#
  eta <- ergm.eta(theta, etamap)
  etagrad <- ergm.etagrad(theta, etamap)
  x <- eta-eta0
  x <- x[!etamap$offsettheta]
  basepred <- xsim %*% x
  prob <- max(basepred)
  prob <- probs*exp(basepred - prob)
  prob <- prob/sum(prob)
  E <- apply(sweep(xsim, 1, prob, "*"), 2, sum)
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
  htmp <- htmp %*% t(etagrad)
  H <- t(htmp) %*% htmp
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
                      varweight=0.5, trustregion, eta0, etamap){
  eta <- ergm.eta(theta, etamap)
  x <- eta-eta0
  basepred <- xsim * x
  maxbase <- max(basepred)
  llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
  llr
}

llik.grad2 <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                       varweight=0.5, trustregion, eta0, etamap){
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
                      varweight=0.5, trustregion, eta0, etamap){ # eqn (5) 
  eta <- ergm.eta(theta, etamap)
  deta <- matrix(eta-eta0,ncol=1) #px1
  basepred <- as.vector(xsim %*% deta) #nx1
  maxbase <- max(basepred)
  sum(xobs * deta) - maxbase - log(sum(probs*exp(basepred-maxbase)))
}

llik.grad3 <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
                       varweight=0.5, eta0, etamap){ #eqn (11)
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
