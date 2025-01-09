#  File tests/testthat/test-miss-dep.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
tolerance <- 4 # Result must be within 4*MCMCSE of truth.

data(florentine)

test_that("Missing data MLE with edge observational constraint", {
  n.active <- network.dyadcount(flomarriage) - network.edgecount(flobusiness)
  e.active <- network.edgecount(flomarriage) - network.edgecount(flomarriage&flobusiness)
  p.active <- e.active/n.active
  theta <- logit(p.active)
  mcmcfit <- suppressWarnings(ergm(flomarriage~edges, constraints = ~fixedas(flobusiness), obs.constraints = ~edges))

  ## Point estimate
  expect_lt(abs(theta-coef(mcmcfit))/sqrt(diag(vcov(mcmcfit, source="estimation"))), tolerance)

  ## Likelihood
  # Relative to NULL model, whose likelihood is defined to be 0.
  llk <- log(p.active)*e.active+log(1-p.active)*(n.active-e.active) - log(.5)*n.active
  expect_lt(abs(llk-logLik(mcmcfit))/sqrt(attr(logLik(mcmcfit), "vcov")), tolerance) 
})
