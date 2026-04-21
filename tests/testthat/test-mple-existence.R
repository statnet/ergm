#  File tests/testthat/test-mple-existence.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

test_that("mple.existence()", {
  pl.nonsep <- list(xmat = matrix(c(1, 1), ncol = 1), zy = c(1, 0))
  pl.sep <- list(xmat = matrix(c(1, -1), ncol = 1), zy = c(1, 0))

  expect_warning(ergm:::mple.existence(pl.nonsep, "lpsolve"), NA)
  expect_warning(ergm:::mple.existence(pl.sep, "lpsolve"), "does not exist")
  expect_warning(ergm:::mple.existence(pl.nonsep, "glpk"), NA)
  expect_warning(ergm:::mple.existence(pl.sep, "glpk"), "does not exist")
})
