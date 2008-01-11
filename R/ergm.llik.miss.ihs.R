#
#   missing data code
#
llik.fun.miss <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
                     penalty=0.5, trustregion=20, eta0, etamap){
  theta.offset <- etamap$theta0
  theta.offset[!etamap$offsettheta] <- theta
  eta <- ergm.eta(theta.offset, etamap)
  x <- eta-eta0
# The next line is right!
# aaa <- sum(xobs * x) - log(sum(probs*exp(xsim %*% x)))
# These lines standardize:
  basepred <- xsim %*% x
  misspred <- xsim.miss %*% x
#
# maxbase <- max(basepred)
# llr <- sum(xobs * x) - maxbase - log(sum(probs*exp(basepred-maxbase)))
#
# alternative based on log-normal approximation
  mb <- sum(basepred*probs)
# vb <- sum(basepred*basepred*probs) - mb*mb
  vb <- sum((basepred-mb)*(basepred-mb)*probs)
  mm <- sum(misspred*probs.miss)
  vm <- sum((misspred-mm)*(misspred-mm)*probs.miss)
# 
# This is the log-likelihood ratio (and not its negative)
#
  llr <- sum(xobs * x) + (mm + penalty*vm) - (mb + penalty*vb)
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
llik.grad.miss <- function(theta, xobs, xsim, probs,  xsim.miss=NULL, probs.miss=NULL,
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
  misspred <- xsim.miss %*% x
  prob.miss <- max(misspred)
  prob.miss <- probs.miss*exp(misspred - prob.miss)
  prob.miss <- prob.miss/sum(prob.miss)
  E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
  llr <- xobs + E.miss-E
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

llik.hessian.miss <- function(theta, xobs, xsim, probs, xsim.miss=NULL, probs.miss=NULL,
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
  misspred <- xsim.miss %*% x
  prob.miss <- max(misspred)
  prob.miss <- probs.miss*exp(misspred - prob.miss)
  prob.miss <- prob.miss/sum(prob.miss)
  E.miss <- apply(sweep(xsim.miss, 1, prob.miss, "*"), 2, sum)
  llr <- xobs + E.miss-E
# 
  htmp <- sweep(sweep(xsim, 2, E, "-"), 1, sqrt(prob), "*")
  htmp <- htmp %*% t(etagrad)
  H <- t(htmp) %*% htmp
  htmp <- sweep(sweep(xsim.miss, 2, E.miss, "-"), 1, sqrt(prob.miss), "*")
  htmp <- htmp %*% etagrad
  H.miss <- t(htmp) %*% htmp
  H <- H.miss-H
  dimnames(H) <- list(namesx, namesx)
  -H
}
