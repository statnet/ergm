#  File tests/testthat/test-scoping.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

test_that("passing components of an ergm() formula into a function", {
  control <- snctrl(MCMLE.maxit = 1, MCMC.burnin = 100, MCMC.interval = 100, seed = 0)

  sim_func <- function(nw, deg, ...) simulate(nw ~ edges + degree(deg), ...)
  fit_func <- function(nw, deg) ergm(nw ~ edges + degree(deg), control = control, eval.loglik=TRUE)
  mynet <- network.initialize(20, directed = FALSE)
  mynet <- sim_func(mynet, 3, coef = c(-1.5, 0))
  ffit <- fit_func(mynet, 3)
  fit <- ergm(mynet ~ edges + degree(3), control = control, eval.loglik=TRUE)

  ffit$call <- ffit$formula <- fit$call <- fit$formula <- NULL
  expect_equal(fit, ffit, ignore_function_env = TRUE, ignore_formula_env = TRUE)
})
