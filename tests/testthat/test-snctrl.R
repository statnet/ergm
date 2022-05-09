#  File tests/testthat/test-snctrl.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
data(florentine)

test_that("snctrl() has at least some of the correct arguments", {
  snctrl_names <- names(formals(snctrl))
  ergm_names <- names(formals(control.ergm))
  simulate.ergm_names <- names(formals(control.simulate.ergm))
  expect_true(all(ergm_names %in% snctrl_names))
  expect_true(all(simulate.ergm_names %in% snctrl_names))
})

test_that("snctrl() works", {
  # Warn on unrecognized argument:
  expect_warning(ctrl <- ergm(flomarriage~edges, control=snctrl(parallel=2, MCMLE.maxit=42, SAN.packagenames="blah", misspelled_argument=7))$control, "^.*misspelled_argument.*$")
  # Set values are correctly assigned in nested lists:
  expect_equal(ctrl$MCMLE.maxit, 42)
  expect_equal(ctrl$SAN$SAN.packagenames, "blah")
  # This includes shared control parameters:
  expect_equal(ctrl$parallel, 2)
  expect_equal(ctrl$SAN$parallel, 2)
})
