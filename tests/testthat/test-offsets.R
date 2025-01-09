#  File tests/testthat/test-offsets.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
o <- options(ergm.eval.loglik=TRUE)

set.seed(0)

data(sampson)
total.theta <- coef(ergm(samplike~edges))
offset.theta <- pi

wald_pval <- function(fit1, fit2, i1=TRUE, i2=TRUE){
  d <- coef(fit1)[i1] - coef(fit2)[i2]
  df <- length(d)
  vcov <- vcov(fit1, sources="estimation")[i1,i1] + vcov(fit2, sources="estimation")[i2,i2]
  pchisq(d %*% solve(vcov) %*% d, df, lower.tail=FALSE)
}

test_that("Linear ERGM with free parameter before offset", {
  e1 <- ergm(samplike~edges+offset(edges), offset.coef=c(pi))
  expect_equal(coef(e1)[1], total.theta-offset.theta, tolerance=0.00001, ignore_attr=TRUE)
  expect_equal(attr(logLik(e1),"df"),1)
})

test_that("Linear ERGM with free parameter after offset", {
  e2 <- ergm(samplike~offset(edges)+edges, offset.coef=c(pi))
  expect_equal(coef(e2)[2], total.theta-offset.theta, tolerance=0.00001, ignore_attr=TRUE)
  expect_equal(attr(logLik(e2),"df"),1)
})

test_that("Linear ERGM with partial offsets", {
  e3 <- ergm(samplike~edges+nodematch("group", diff=TRUE, levels=-2))
  e3a <- ergm(samplike~edges+offset(nodematch("group", diff=TRUE), 2), offset.coef=0)
  e3b <- ergm(samplike~edges+offset(nodematch("group", diff=TRUE), c(FALSE,TRUE,FALSE)), offset.coef=0)
  expect_equal(coef(e3a)[-3], coef(e3), tolerance=0.00001)
  expect_equal(coef(e3b)[-3], coef(e3), tolerance=0.00001)
  expect_equal(logLik(e3a), logLik(e3), tolerance=0.00001)
  expect_equal(logLik(e3b), logLik(e3), tolerance=0.00001)
})

test_that("Curved ERGM with partial offsets", {
  e4 <- ergm(samplike~edges+gwesp(0.25, fix=TRUE), control=control.ergm(seed=0,MCMLE.maxit=2))
  e4a <- ergm(samplike~edges+offset(gwesp(),c(FALSE,TRUE)), offset.coef=0.25, control=control.ergm(seed=0,MCMLE.maxit=2))
  expect_gt(wald_pval(e4a, e4, -3, TRUE), 0.01)
  expect_equal(logLik(e4a), logLik(e4), tolerance=0.01, ignore_attr=TRUE)
})


options(o)
