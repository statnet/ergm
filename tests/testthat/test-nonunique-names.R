#  File tests/testthat/test-nonunique-names.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
data(samplk)
samplk2 %e% "a" <- 1
samplk3 %e% "a" <- 1
test_that("MCMC diagnostics produced even when names are not unique", {
  fit <- ergm(samplk1~edgecov(samplk2,"a")+edgecov(samplk3,"a"), control=control.ergm(force.main=TRUE, MCMLE.maxit=1, MCMC.burnin=1, MCMC.interval=1), eval.loglik=FALSE)
  mcmc.diagnostics(fit) |> expect_no_condition()
})
