context("test-nonunique-names.R")
data(samplk)
samplk2 %e% "a" <- 1
samplk3 %e% "a" <- 1
test_that("MCMC diagnostics produced even when names are not unique", {
  fit <- ergm(samplk1~edgecov(samplk2,"a")+edgecov(samplk3,"a"), control=control.ergm(force.main=TRUE, MCMLE.maxit=1, MCMC.burnin=1, MCMC.interval=1), eval.loglik=FALSE)
  mcmc.diagnostics(fit)
})
