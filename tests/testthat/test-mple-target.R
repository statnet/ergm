context("test-mple-target.R")
n<-500
base.net <- network.initialize(n=n,directed=FALSE)
norm.stats<-c(.7,.1,.5)
target.stats<-norm.stats*n
print(target.stats)
cat("Structural check:\n")
cat("Mean degree:", norm.stats[1]*2,".\n")
cat("Average degree among nodes with degree 2 or higher:", (2*norm.stats[1]-norm.stats[3])/(1-norm.stats[2]-norm.stats[3]),".\n")

ergm.ts.fit<-ergm(base.net~edges+degree(c(0,1)),target.stats=n*norm.stats,estimate="MPLE")

test_that("internal SAN call achieves the target statistics", {
  expect_equivalent(summary(ergm.ts.fit$newnetwork~edges+degree(c(0,1)),estimate="MPLE"),target.stats)
})

test_that("estimate with target.stats matches that with LHS", {
  ergm.fit<-ergm(ergm.ts.fit$network~edges+degree(c(0,1)),estimate="MPLE")
  expect_equal(coef(ergm.fit),coef(ergm.ts.fit))
})

test_that("simulating from the MPLE target statistics fit", {
  ergm.sim<-simulate(ergm.ts.fit,nsim=10,output="stats", control=control.simulate.ergm(MCMC.burnin=10,MCMC.interval=1))
})
