#  File tests/testthat/test-mple-existence.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

test_that("MPLE.check control accepts solver and skip values", {
  expect_equal(ergm:::resolve_MPLE_check(TRUE), "glpk")
  expect_equal(ergm:::resolve_MPLE_check("glpk"), "glpk")
  expect_equal(ergm:::resolve_MPLE_check("lpsolve"), "lpsolve")
  expect_equal(ergm:::resolve_MPLE_check(FALSE), "skip")
  expect_equal(ergm:::resolve_MPLE_check("skip"), "skip")
  expect_error(ergm:::resolve_MPLE_check("bogus"), "should be one of")
})

test_that("mple.existence() supports lpsolve solver selection", {
  pl.nonsep <- list(xmat=matrix(c(1, 1), ncol=1), zy=c(1, 0))
  pl.sep <- list(xmat=matrix(c(1, -1), ncol=1), zy=c(1, 0))

  expect_warning(ergm:::mple.existence(pl.nonsep, solver="lpsolve"), NA)
  expect_warning(ergm:::mple.existence(pl.sep, solver="lpsolve"), "does not exist")
})
