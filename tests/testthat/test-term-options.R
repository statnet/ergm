#  File tests/testthat/test-term-options.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
context("test-term-options.R")

data(florentine)
old.opts1 <- options(ergm.eval.loglik=TRUE)
times <- 2

test_that("summary() of a term that takes options",{
  expect_equivalent(summary(flomarriage~.edges_times, term.options=list(times=2)),
                    summary(flomarriage~edges)*times)
})

e1 <- ergm(flomarriage~.edges_times, control=control.ergm(term.options=list(times=2)))
test_that("ergm() MPLE of a term that takes options",{
  expect_equivalent(coef(e1),coef(ergm(flomarriage~edges))/times)
})

e2 <- ergm(flomarriage~.edges_times, control=control.ergm(force.main=TRUE,seed=0,term.options=list(times=2)))
test_that("ergm() MCMLE of a term that takes options",{
  expect_equivalent(coef(e2),coef(e1), tolerance=.005)
  expect_equivalent(logLik(e2),logLik(e1), tolerance=.005)
})

test_that("gof() of a model that had term options",{
  gof(e2, control=control.gof.ergm(nsim=10))
})

options(ergm.eval.loglik=FALSE)
e3 <- ergm(flobusiness~.edges_times, target.stats=as.vector(summary(flomarriage~edges)/2), control=control.ergm(seed=0,force.main=TRUE,term.options=list(times=2)))
test_that("ergm() for a model with term options and mean-value parameters and log-likelihood calculation off via a global option",{
  expect_equivalent(coef(e2),coef(e1), tolerance=.005)
  expect_error(logLik(e3), "Log-likelihood was not estimated for this fit.")
})

old.opts2 <- options(ergm.term=list(times=2))
test_that("summary() with term options globally set",{
  expect_equivalent(summary(samplike~.edges_times),summary(samplike~edges)*times)
})


old.opts3 <- options(ergm.term=list(times=1))
test_that("control parameter overrides global setting",{
  expect_equivalent(summary(samplike~.edges_times, term.options=list(times=2)),summary(samplike~edges)*times)
})
options(old.opts3)
options(old.opts2)
options(old.opts1)
