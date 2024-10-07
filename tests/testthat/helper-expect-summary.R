#  File tests/testthat/helper-expect-summary.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
expect_summary <- function(s, e, value, coefficients, tolerance=0.001) {
  expect_equal(s, value, tolerance=tolerance, ignore_attr=is.null(names(value)))
  expect_equal(coef(e)[1:length(coefficients)], coefficients, tolerance=tolerance, ignore_attr=is.null(names(coefficients)))
}
