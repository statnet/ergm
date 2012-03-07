#  File ergm/tests/miss_tests.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)

theta0err<--1 # Perturbation in the initial values
maxit<-20 # Maximum number of iterations
tolerance<-0.01 # Result must be within 1% of truth.

n<-20 # Number of nodes
b<-3 # Bipartite split

d<-.1 # Density
m<-.1 # Missingness rate

logit<-function(p) log(p/(1-p))

cat("n=",n,", density=",d,", missing=",m,"\n",sep="")
mk.missnet<-function(n,d,m,directed=TRUE,bipartite=0){
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
  
  mplefit<-ergm(y~edges)
  mpleOK<-all.equal(truth, coef(mplefit), check.attributes=FALSE, tolerance=tolerance)
  cat("MPLE estimate =", coef(mplefit), if(isTRUE(mpleOK)) "OK" else mpleOK,"\n")

  mcmcfit<-ergm(y~edges, control=control.ergm(force.main=TRUE, MCMC.interval=ceiling(n^(3/2)), MCMLE.maxit=maxit, 
               init=truth+theta0err))
  mcmcOK<-all.equal(truth, coef(mcmcfit), check.attributes=FALSE, tolerance=tolerance)
  cat("MCMCMLE estimate =", coef(mcmcfit), if(isTRUE(mcmcOK)) "OK" else mcmcOK,"\n")
  
  return(isTRUE(mpleOK) && isTRUE(mcmcOK))
}

# Directed
cat("\n\nDirected Network\n")
set.seed(123)
y<-mk.missnet(n, d, m, TRUE, 0)
stopifnot(run.miss.test(y))

# Undirected
cat("\n\nUndirected Network\n")
set.seed(456)
y<-mk.missnet(n, d, m, FALSE, 0)
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
y.miss <- simulate(y~edges, coef=logit(0.1))
y[as.edgelist(y.miss)] <- NA

cat("Network statistics:\n")
print(summary(y~edges+gwesp(0.5)))
truth<-correct.edges.theta(y)
cat("Correct estimate =",truth,"\n")

set.seed(654)
mcmcfit<-ergm(y~edges+gwesp(0.5), control=control.ergm(MCMLE.maxit=5))
summary(mcmcfit)
stopifnot(abs(coef(mcmcfit)[1]-truth)/sqrt(mcmcfit$covar[1])<2)
