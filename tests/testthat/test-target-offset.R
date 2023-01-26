#  File tests/testthat/test-target-offset.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
data(florentine)

test_that("target+offset in a non-curved ERGM", {
  expect_warning(fit <- ergm(flomarriage~offset(edges)+edges+degree(1)+offset(degree(0)),target.stats=summary(flomarriage~edges+degree(1)),
                             offset.coef=c(0,-0.25),control=control.ergm(init=c(0,-1.47,0.462,-0.25),MPLE.nonident.tol=0,MCMLE.nonident.tol=0), eval.loglik=TRUE),
                 "^Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood.*")
  expect_warning(expect_error(summary(fit),NA),NA)
  expect_warning(expect_error(mcmc.diagnostics(fit),NA),NA)
})

test_that("target+offset in a curved ERGM", {
  set.seed(10)
  suppressWarnings(expect_warning(
      fit <- ergm(flomarriage~offset(edges)+edges+gwdegree()+degree(0)+offset(degree(1)),target.stats=summary(flomarriage~edges+gwdegree(fix=FALSE)+degree(0)), offset.coef=c(0,-0.25), control=control.ergm(MCMLE.termination="none",MCMLE.maxit=3,MPLE.nonident.tol=0,MCMLE.nonident.tol=0), eval.loglik=TRUE),
      "^Using target.stats for a model with offset terms may produce an inaccurate estimate of the log-likelihood.*"))
  expect_warning(expect_error(summary(fit),NA),NA)
  expect_warning(expect_error(mcmc.diagnostics(fit),NA),NA)
})
