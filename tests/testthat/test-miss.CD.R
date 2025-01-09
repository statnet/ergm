#  File tests/testthat/test-miss.CD.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

attach(MLE.tools)

theta0err<--1 # Perturbation in the initial values
maxit<-60 # Maximum number of iterations
tolerance<-0.01 # Result must be within 1% of truth.
tolerance.CD<-0.15 # Result must be within 15% of truth.

n<-20 # Number of nodes
b<-7 # Bipartite split

d<-.1 # Density
m<-.1 # Missingness rate

cat("n=",n,", density=",d,", missing=",m,"\n",sep="")

run.miss.test<-function(y){
  truth<-edges.theta(y)

  ### Needs more work.
  ## cdfit<-ergm(y~edges, estimate="CD")
  ## cdOK<-all.equal(truth, coef(cdfit), ignore_attr=TRUE, tolerance=tolerance.CD)
  ## cat("CD estimate =", coef(cdfit), if(isTRUE(cdOK)) "OK" else cdOK,"\n")

  cd2fit<-ergm(y~edges, control=control.ergm(CD.nsteps=50, MCMC.samplesize=100), estimate="CD")
  expect_equal(truth, coef(cd2fit), ignore_attr=TRUE, tolerance=tolerance.CD)
}

# Directed
test_that("directed network", {
  set.seed(123)
  y<-mk.missnet(n, d, m, TRUE, FALSE)
  run.miss.test(y)
})

# Undirected
test_that("undirected network", {
  set.seed(456)
  y<-mk.missnet(n, d, m, FALSE, FALSE)
  run.miss.test(y)
})

# Bipartite Undirected
test_that("bipartite undirected network", {
  set.seed(789)
  y<-mk.missnet(n, d, m, FALSE, b)
  run.miss.test(y)
})

# Add the curved+missing test here for now
test_that("curved+missing", {
  set.seed(321)
  n <- 30
  y <- network.initialize(n, directed=FALSE) # Create an empty network
  y <- simulate(y~edges, coef=logit(0.12), control=control.simulate(MCMC.burnin=2*n^2))
  y.miss <- simulate(y~edges, coef=logit(0.01))
  y[as.edgelist(y.miss)] <- NA

  cat("Network statistics:\n")
  print(summary(y~edges+gwesp()))
  truth<-edges.theta(y)
  cat("Correct estimate =",truth,"\n")

  set.seed(654)
  cdfit<-ergm(y~edges+gwesp(), estimate="CD", control=control.ergm(CD.nsteps=50, MCMC.samplesize=100))
  summary(cdfit)
  expect_lt(abs(coef(cdfit)[1]-truth)/sqrt(cdfit$covar[1]), 2)
})

detach(MLE.tools)
