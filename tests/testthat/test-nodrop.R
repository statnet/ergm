#  File tests/testthat/test-nodrop.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
### Tests to make sure drop=FALSE works.

library(statnet.common)
opttest({

data(sampson)
test_that("shouldn't need to drop", {
  # MPLE
  summary(ergm(samplike~edges, control=control.ergm(drop=FALSE)))
  # MCMC
  summary(ergm(samplike~edges, control=control.ergm(drop=FALSE, force.main=TRUE, MCMLE.maxit=3)))
})

test_that("empty network", {
  y0 <- network.initialize(10)
  # MPLE
  summary(ergm(y0~edges, control=control.ergm(drop=FALSE)))
  # MCMC
  summary(ergm(y0~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=3)))
})

test_that("full network", {
  y1 <- as.network(matrix(1,10,10))
  # MPLE
  summary(ergm(y1~edges, control=control.ergm(drop=FALSE)))
  # MCMC
  summary(ergm(y1~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=3)))
})

}, "drop disabled")
