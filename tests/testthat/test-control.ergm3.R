#  File tests/testthat/test-control.ergm3.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
test_that("control.ergm3() defaults differ", {
  expect_equal(lapply(diff(control.ergm(), control.ergm3())[c("MCMLE.termination", "MCMLE.effectiveSize")], `[[`, "y"),
               list(MCMLE.termination = "Hummel", MCMLE.effectiveSize = NULL))
})

test_that("control.ergm3() is otherwise identical", {
  ## Two distinct tests are needed because the first one is a more
  ## natural use, but it affects the environment in which some
  ## parameters that are functions are defined.

  expect_equal(control.ergm(), control.ergm3(MCMLE.termination = "confidence", MCMLE.effectiveSize = 64), ignore_function_env = TRUE)

  expect_equal(control.ergm(), control.ergm3(MCMLE.termination = c("confidence", "Hummel", "Hotelling", "precision", "none"), MCMLE.effectiveSize = 64))
})
