#  File tests/valued_sim.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
load("testnet3d.RData")

set.seed(0)

## DiscUnif-reference
cat("======== Discrete-uniform-reference ERGM with minimum of -1 and maxium of 5\n")

f <- function(x, q, a, b) exp(q*x)/sum(exp(q*(a:b)))

E <- function(q, a, b) sum((a:b)*exp(q*(a:b))/sum(exp(q*(a:b))))

V <- function(q, a, b) sum(((a:b)-E(q,a,b))^2*exp(q*(a:b))/sum(exp(q*(a:b))))

a <- -1; b <- 5; coefs <- c(runif(1,-3,3),0)

for(coef in coefs){
  cat("==== statsonly=TRUE, coef=",coef,"\n",sep="")
  s <- simulate(testnet3d~sum, nsim=1000, reference=~DiscUnif(a,b), response="w", coef=coef, statsonly=TRUE, control=control.simulate(MCMC.burnin=100))
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}

  
  cat("==== statsonly=FALSE, coef=",coef,"\n",sep="")
  s.full<-simulate(testnet3d~sum, nsim=1000, reference=~DiscUnif(a,b), response="w", coef=coef, statsonly=FALSE, control=control.simulate(MCMC.burnin=100))
  s <- sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w")),simplify=TRUE)
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}
}
  
## Unif-reference
cat("======== Continuous-uniform-reference ERGM with minimum of -1 and maxium of 5\n")

E <- function(q, a, b) if(isTRUE(all.equal(q,0))) (b+a)/2 else ((b*q-1)*exp(b*q)-(a*q-1)*exp(a*q))/(q*exp(b*q)-q*exp(a*q))
V <- function(q, a, b) if(isTRUE(all.equal(q,0))) (b-a)^2/12 else (((-b^2+2*a*b-a^2)*q^2-2)*exp(b*q+a*q)+exp(2*b*q)+exp(2*a*q))/(-2*q^2*exp(b*q+a*q)+q^2*exp(2*b*q)+q^2*exp(2*a*q))

a <- -1; b <- 5; coefs <- c(runif(1,-3,3),0)

for(coef in coefs){
  cat("==== statsonly=TRUE, coef=",coef,"\n",sep="")
  s <- simulate(testnet3d~sum, nsim=1000, reference=~Unif(a,b), response="w", coef=coef, statsonly=TRUE, control=control.simulate(MCMC.burnin=100))
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}

  
  cat("==== statsonly=FALSE, coef=",coef,"\n",sep="")
  s.full<-simulate(testnet3d~sum, nsim=1000, reference=~Unif(a,b), response="w", coef=coef, statsonly=FALSE, control=control.simulate(MCMC.burnin=100))
  s <- sapply(s.full, function(x) sum(as.matrix(x,m="a",a="w")),simplify=TRUE)
  test <- approx.hotelling.diff.test(s/6,mu0=E(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (mean).")}
  test <- approx.hotelling.diff.test((s-E(coef,a,b)*6)^2/6,mu0=V(coef,a,b))
  if(test$p.value<0.001) {print(test); stop("Simulation test failed (variance).")}
}
}, "continuous uniform reference")
