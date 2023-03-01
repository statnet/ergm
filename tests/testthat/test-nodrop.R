#  File tests/testthat/test-nodrop.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
### Tests to make sure drop=FALSE works.

o <- options(ergm.eval.loglik=FALSE)

data(sampson)
test_that("shouldn't need to drop", {
  # MPLE
  expect_warning(ergm(samplike~edges, control=control.ergm(drop=FALSE)), NA)
  # MCMC
  expect_warning(ergm(samplike~edges, control=control.ergm(drop=FALSE, force.main=TRUE, MCMLE.maxit=2)), NA)
})

test_that("empty network", {
  y0 <- network.initialize(10)
  # MPLE
  expect_warning(expect_warning(expect_warning(
    ergm(y0~edges, control=control.ergm(drop=FALSE)),
    "The MPLE does not exist!"),
    "Network is empty and no target stats are specified."),
    "Observed statistic.s. edges are at their smallest attainable values and drop=FALSE. The MLE is poorly defined.")
  # MCMC
  expect_warning(expect_warning(
    ergm(y0~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=2)),
    "Network is empty and no target stats are specified."),
    "Observed statistic.s. edges are at their smallest attainable values and drop=FALSE. The MLE is poorly defined.")
})

test_that("full network", {
  y1 <- as.network(matrix(1,10,10))
  # MPLE
  expect_warning(expect_warning(
    ergm(y1~edges, control=control.ergm(drop=FALSE)),
    "The MPLE does not exist!"),
    "Observed statistic.s. edges are at their greatest attainable values and drop=FALSE. The MLE is poorly defined.")

  # MCMC
  expect_warning(
    ergm(y1~edges, control=control.ergm(drop=FALSE, force.main=TRUE, init=0, MCMLE.maxit=2)),
    "Observed statistic.s. edges are at their greatest attainable values and drop=FALSE. The MLE is poorly defined.")
})

options(o)
