#  File tests/testthat/test-shrink-into-CH.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
test_that("shrink_into_CH() works in both Rglpk and lpSolveAPI modes and produces identical results", {
  set.seed(0)
  p <- matrix(rnorm(200), 20)
  M <- matrix(rnorm(2000), 200)
  m <- colMeans(M)
  expect_equal(shrink_into_CH(p, M, m, solver="glpk"), shrink_into_CH(p, M, m, solver="lpsolve"))
})
