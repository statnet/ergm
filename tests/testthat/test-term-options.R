#  File tests/testthat/test-term-options.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

data(florentine)
old.opts1 <- options(ergm.eval.loglik=TRUE)
times <- 2

test_that("summary() of a term that takes options",{
  expect_equal(summary(flomarriage~.edges_times, term.options=list(times=2)),
               summary(flomarriage~edges)*times, ignore_attr=TRUE)
})

e1 <- ergm(flomarriage~.edges_times, control=control.ergm(term.options=list(times=2)))
test_that("ergm() MPLE of a term that takes options",{
  expect_equal(coef(e1),coef(ergm(flomarriage~edges))/times, ignore_attr=TRUE)
})

e2 <- ergm(flomarriage~.edges_times, control=control.ergm(force.main=TRUE,seed=0,term.options=list(times=2)))
test_that("ergm() MCMLE of a term that takes options",{
  expect_equal(coef(e2),coef(e1), tolerance=.05, ignore_attr=TRUE)
  expect_equal(logLik(e2),logLik(e1), tolerance=.05, ignore_attr=TRUE)
})

test_that("gof() of a model that had term options",{
  expect_error(gof(e2, control=control.gof.ergm(nsim=10)), NA)
})

old.opts2 <- options(ergm.eval.loglik=FALSE)
e3 <- ergm(flobusiness~.edges_times, target.stats=as.vector(summary(flomarriage~edges)/2), control=control.ergm(seed=0,force.main=TRUE,term.options=list(times=2)))
test_that("ergm() for a model with term options and mean-value parameters and log-likelihood calculation off via a global option",{
  expect_equal(coef(e2),coef(e1), tolerance=.05, ignore_attr=TRUE)
  expect_error(logLik(e3), "Log-likelihood was not estimated for this fit.")
})

data(sampson)
old.opts3 <- options(ergm.term=list(times=2))
test_that("summary() with term options globally set",{
  expect_equal(summary(samplike~.edges_times),summary(samplike~edges)*times, ignore_attr=TRUE)
})


old.opts4 <- options(ergm.term=list(times=1))
test_that("control parameter overrides global setting",{
  expect_equal(summary(samplike~.edges_times, term.options=list(times=2)),summary(samplike~edges)*times, ignore_attr=TRUE)
})

options(old.opts4)
options(old.opts3)
options(old.opts2)
options(old.opts1)
