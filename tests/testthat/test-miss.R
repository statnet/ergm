#  File tests/testthat/test-miss.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

attach(MLE.tools)

theta0err<- 1 # Perturbation in the initial values
tolerance <- 4 # Result must be within 4*MCMCSE of truth.
bridge.target.se <- 0.005 # Log-likelihood MCMC standard error must be within this.

n<-20 # Number of nodes
b<-3 # Bipartite split

d<-.1 # Density
m<-.05 # Missingness rate

cat("n=",n,", density=",d,", missing=",m,"\n",sep="")

run.miss.test<-function(y){
  theta <- edges.theta(y)
  cat("Correct estimate =",theta,"with log-likelihood",edges.llk(y),".\n")

  mplefit<-ergm(y~edges, eval.loglik=TRUE)
  mple.theta.OK<-all.equal(theta,coef(mplefit),check.attributes=FALSE)
  mple.llk.OK<-all.equal(edges.llk(y, coef(mplefit)),
                         as.vector(logLik(mplefit)),check.attributes=FALSE)
  cat("MPLE estimate =", coef(mplefit),"with log-likelihood",logLik(mplefit), if(isTRUE(mple.theta.OK)&&isTRUE(mple.llk.OK)) "OK.","\n")

  mcmcfit<-ergm(y~edges, control=snctrl(force.main=TRUE, init=theta+theta0err, bridge.target.se=bridge.target.se), verbose=TRUE, eval.loglik=TRUE)
  mcmc.diagnostics(mcmcfit)
  mcmc.theta.OK<-abs(theta-coef(mcmcfit))/sqrt(diag(vcov(mcmcfit, source="estimation")))
  mcmc.llk.OK<-abs(edges.llk(y, coef(mcmcfit))-logLik(mcmcfit))/bridge.target.se

  cat("MCMCMLE estimate =", coef(mcmcfit),"with log-likelihood",logLik(mcmcfit), if(mcmc.theta.OK<tolerance&&mcmc.llk.OK<tolerance) "OK.","\n")

  expect_true(is.na(mplefit))
  expect_true(is.na(mcmcfit))
  expect_true(anyNA(mplefit))
  expect_true(anyNA(mcmcfit))
  expect_true(is.dyad.independent(mplefit))
  expect_true(is.dyad.independent(mcmcfit))
  expect_true(isTRUE(mple.theta.OK))
  expect_lt(mcmc.theta.OK, tolerance)
  expect_true(isTRUE(mple.llk.OK))
  expect_lt(mcmc.llk.OK, tolerance)
}

# Directed
test_that("directed network", {
  set.seed(123)
  y<-mk.missnet(n, d, m, TRUE, FALSE)
  run.miss.test(y)
})

# Undirected
test_that("undirected Network", {
  set.seed(456)
  y<-mk.missnet(n, d, m, FALSE, FALSE)
  run.miss.test(y)
})

# Bipartite Undirected
test_that("Bipartite Undirected Network", {
  set.seed(0)
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
  mcmcfit<-ergm(y~edges+gwesp(), control=control.ergm(MCMLE.maxit=5))
  summary(mcmcfit)
  expect_lt(abs(coef(mcmcfit)[1]-truth)/sqrt(mcmcfit$covar[1]), tolerance)
})

detach(MLE.tools)
