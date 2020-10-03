#  File tests/miss_tests.CD.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
theta0err<--1 # Perturbation in the initial values
maxit<-60 # Maximum number of iterations
tolerance<-0.01 # Result must be within 1% of truth.
tolerance.CD<-0.1 # Result must be within 10% of truth.

n<-20 # Number of nodes
b<-7 # Bipartite split

d<-.1 # Density
m<-.1 # Missingness rate

logit<-function(p) log(p/(1-p))

cat("n=",n,", density=",d,", missing=",m,"\n",sep="")
mk.missnet<-function(n,d,m,directed=TRUE,bipartite=FALSE){
  y<-network.initialize(n, directed=directed, bipartite=bipartite)
  y<-simulate(y~edges, coef=logit(d), control=control.simulate(MCMC.burnin=2*n^2))
  if(m>0){
    y.miss<-simulate(y~edges, coef=logit(m))
    y[as.edgelist(y.miss)]<-NA
  }
  y
}

correct.edges.theta<-function(y){
  e<-summary(y~edges)
  d<-network.dyadcount(y)
  m<-network.naedgecount(y)

  logit(e/d)
}


run.miss.test<-function(y){
  truth<-correct.edges.theta(y)
  cat("Correct estimate =",truth,"\n")

  ### Needs more work.
  ## cdfit<-ergm(y~edges, estimate="CD")
  ## cdOK<-all.equal(truth, coef(cdfit), check.attributes=FALSE, tolerance=tolerance.CD)
  ## cat("CD estimate =", coef(cdfit), if(isTRUE(cdOK)) "OK" else cdOK,"\n")

  cd2fit<-ergm(y~edges, control=control.ergm(CD.nsteps=50, MCMC.samplesize=100), estimate="CD")
  cd2OK<-all.equal(truth, coef(cd2fit), check.attributes=FALSE, tolerance=tolerance.CD)
  cat("CD2 estimate =", coef(cd2fit), if(isTRUE(cd2OK)) "OK" else cd2OK,"\n")

  
  return(
      ## isTRUE(cdOK) &&                   
      isTRUE(cd2OK))
}

# Directed
cat("\n\nDirected Network\n")
set.seed(123)
y<-mk.missnet(n, d, m, TRUE, FALSE)
stopifnot(run.miss.test(y))

# Undirected
cat("\n\nUndirected Network\n")
set.seed(456)
y<-mk.missnet(n, d, m, FALSE, FALSE)
stopifnot(run.miss.test(y))  

# Bipartite Undirected
cat("\n\nBipartite Undirected Network\n")
set.seed(789)
y<-mk.missnet(n, d, m, FALSE, b)
stopifnot(run.miss.test(y))

# Add the curved+missing test here for now

set.seed(321)
n <- 50
y <- network.initialize(n, directed=FALSE) # Create an empty network
y <- simulate(y~edges, coef=logit(0.12), control=control.simulate(MCMC.burnin=2*n^2))
y.miss <- simulate(y~edges, coef=logit(0.01))
y[as.edgelist(y.miss)] <- NA

cat("Network statistics:\n")
print(summary(y~edges+gwesp(0.5)))
truth<-correct.edges.theta(y)
cat("Correct estimate =",truth,"\n")

set.seed(654)
cdfit<-ergm(y~edges+gwesp(0.5), estimate="CD", control=control.ergm(CD.nsteps=50, MCMC.samplesize=100))
summary(cdfit)
stopifnot(abs(coef(cdfit)[1]-truth)/sqrt(cdfit$covar[1])<2)
}, "CD missing data")
