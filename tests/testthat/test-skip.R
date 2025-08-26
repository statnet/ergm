#  File tests/testthat/test-skip.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

test_that("simulating with skipping works", {
  nw <- network.initialize(5, dir = FALSE)
  nw[cbind(1:4, 2:5)] <- 1
  coef <- coef(ergm(nw ~ edges))

  stats <- simulate(nw ~ edges + degree(1), monitor = ~degree(0:4),
                    coef = c(coef, NaN), output = "stats", nsim = 1000,
                    control = snctrl(MCMC.burnin = 1, MCMC.interval = 1))

  expect_equal(apply(stats[, 3:7], 2, range) != 0,
               c(FALSE, TRUE,
                 FALSE, FALSE,
                 FALSE, TRUE,
                 FALSE, TRUE,
                 FALSE, TRUE), ignore_attr = TRUE)
})


test_that("estimation with skipping works", {
  nw <- network.initialize(5, dir = FALSE)
  nw %v% "a" <- rep(c(FALSE, TRUE), 2:3)
  nw[3, 1:5] <- 1

  truth <- c(edges = logit(4 / 9), `offset(nodematch.a)` = NaN)

  MPLE.glm.fit <- ergm(nw ~ edges + offset(nodematch("a", levels = 1)), offset.coef = NaN)
  expect_equal(coef(MPLE.glm.fit), truth)

  MPLE.logitreg.fit <- ergm(nw ~ edges + offset(nodematch("a", levels = 1)), offset.coef = NaN,
                            control = snctrl(MPLE.type = "logitreg"))
  expect_equal(coef(MPLE.logitreg.fit), truth)

  MCMLE.fit <- ergm(nw ~ edges + offset(nodematch("a", levels = 1)), offset.coef = NaN,
                    control = snctrl(init.method = "zeros", force.main = TRUE),
                    eval.loglik = FALSE)
  expect_within_mc_err(MCMLE.fit, truth[1], 1L)
})
