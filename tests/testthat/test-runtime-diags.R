#  File tests/testthat/test-runtime-diags.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

o <- options(ergm.eval.loglik=FALSE)
data(florentine)

test_that("runtime diagnostics", {
  expect_error(ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
               control=control.ergm(MCMC.runtime.traceplot=TRUE)), NA)
})
options(o)
