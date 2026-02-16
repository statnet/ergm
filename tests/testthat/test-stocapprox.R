#  File tests/testthat/test-stocapprox.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

options(ergm.eval.loglik=FALSE)

unloadNamespace("ergm.count")

data(florentine)

test_that("Stochastic Approximation produces similar results to MCMLE (linear ERGM)",{
  set.seed(2)

  mod.sa <- ergm(flomarriage~edges+triangle,control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+triangle)

  expect_equal(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1, ignore_attr=TRUE)
})

test_that("Stochastic Approximation produces similar results to MCMLE (curved ERGM)",{
  set.seed(3)

  mod.sa <- ergm(flomarriage~edges+gwdegree(),control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+gwdegree())

  expect_equal(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1, ignore_attr=TRUE)
})


test_that("Stochastic Approximation produces similar results to MCMLE (valued ERGM)",{
  set.seed(2)

  flomarriage %e% "w" <- 1

  mod.sa <- ergm(flomarriage~edges+transitiveweights, response="w", reference=~DiscUnif(0,1), control=control.ergm(main.method="Stochastic-Approximation"))

  mod.mcmle <- ergm(flomarriage~edges+transitiveweights, response="w", reference=~DiscUnif(0,1))

  expect_equal(coef(mod.sa), coef(mod.mcmle), tolerance = 0.1, ignore_attr=TRUE)
})

library(ergm.count)
