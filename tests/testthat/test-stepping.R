#  File tests/testthat/test-stepping.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

library(statnet.common)
opttest({

####Load the data (provided in the package):
data(ecoli, package="statnet.data")
form <- ecoli2 ~ edges + degree(2:5) + gwdegree(0.25, fixed = T)

test_that("stepping test", {
  m2<-ergm(formula=form, verbose=FALSE,
          control=control.ergm(main.method="Stepping", Step.samplesize=100, Step.gridsize=10000,
          MCMLE.metric="lognormal", MCMC.samplesize=1000, MCMC.burnin=1e+4, MCMC.interval=1000,
          seed=12345))
  expect_gte(m2$iterations, 5)
  expect_lte(m2$iterations, 25)
})

}, "Stepping test")
