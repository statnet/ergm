#  File tests/testthat/test-mple-offset.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

test_that("MPLE + offset", {
  set.seed(0)

  options(ergm.eval.loglik=FALSE)
  data(florentine)
  boo<-flomarriage
  boo[1:3,]<-0
  foo <- suppressWarnings(
    ergm(flomarriage~edges+offset(edgecov(boo))+gwesp(0.25,fixed=T),offset.coef=20))
  expect_lte(max(abs(coef(foo)), na.rm=T), 20)
})
