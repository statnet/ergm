#  File tests/testthat/test-term-errors.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
test_that("ergm_Init_abort() can locate the term from which it had been called.", {
  data(faux.mesa.high)
  expect_error(summary(faux.mesa.high~nodefactor("Blah")), ".*term.*\\bnodefactor\\b.*in package.*\\bergm\\b.*")
})
