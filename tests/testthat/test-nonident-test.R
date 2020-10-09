#  File tests/testthat/test-nonident-test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
o <- options(ergm.eval.loglik=FALSE)

data(florentine)

test_that("Nonidentifiable model produces a warning.", {
  expect_warning(ergm(flomarriage~edges+nodecov(~wealth)+nodecov(~-wealth/2+1)),".*nodecov\\.-wealth/2\\+1.*\\bnonidentifiable\\b.*")
  expect_warning(ergm(flomarriage~edges+nodecov(~wealth)+nodecov(~-wealth/2+1), control=control.ergm(init.method="CD")),".*nodecov\\.-wealth/2\\+1.*\\bnonidentifiable\\b.*")
  expect_warning(ergm(flomarriage~edges+nodecov(~wealth)+nodecov(~-wealth/2+1)+gwesp(.1, fixed=FALSE), control=control.ergm(MCMLE.maxit=1)),".*nodecov\\.-wealth/2\\+1.*\\bnonidentifiable\\b.*")
})

test_that("Model identifiable only due to offsets does not.", {
  expect_warning(ergm(flomarriage~edges+offset(nodecov(~wealth))+nodecov(~-wealth/2+1), offset.coef=-1),NA)
  expect_warning(ergm(flomarriage~offset(edges)+nodecov(~wealth)+nodecov(~-wealth/2+1), offset.coef=-1),NA)
  expect_warning(ergm(flomarriage~offset(edges)+nodecov(~wealth)+nodecov(~-wealth/2+1), offset.coef=-1, control=control.ergm(init.method="CD")),NA)
  expect_warning(ergm(flomarriage~edges+offset(nodecov(~wealth))+nodecov(~-wealth/2+1)+gwesp(.1, fixed=FALSE), offset.coef=-1, control=control.ergm(MCMLE.maxit=1)),NA)
})

options(o)
