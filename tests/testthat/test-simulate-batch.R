#  File tests/testthat/test-simulate-batch.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
test_that("simulate.formula() returns the same number of networks and stats regardless of batch size and they are (on average) similar", {
  nw0 <- network.initialize(10)

  set.seed(1)

  unbatched <- simulate(nw0 ~ edges, nsim = 100, coef = -1, control = snctrl(), simplify = FALSE)
  batched1 <- simulate(nw0 ~ edges, nsim = 100, coef = -1, control = snctrl(MCMC.batch = 1), simplify = FALSE)
  batched2 <- simulate(nw0 ~ edges, nsim = 100, coef = -1, control = snctrl(MCMC.batch = 2), simplify = FALSE)

  expect_equal(lapply(attr(unbatched, "stats"), dim), lapply(attr(batched1, "stats"), dim))
  expect_equal(lapply(attr(unbatched, "stats"), dim), lapply(attr(batched2, "stats"), dim))

  expect_equal(lengths(unbatched), lengths(batched1))
  expect_equal(lengths(unbatched), lengths(batched2))

  expect_gte(suppressWarnings(ks.test(as.matrix(attr(unbatched, "stats")), as.matrix(attr(batched1, "stats")))$p.value), .05)
  expect_gte(suppressWarnings(ks.test(as.matrix(attr(unbatched, "stats")), as.matrix(attr(batched2, "stats")))$p.value), .05)
})

test_that("simulate.formula() returns the same number of networks and stats regardless of batch size and they are (on average) similar, with parallel", {
  nw0 <- network.initialize(10)

  # NB: setup_strategy= is a workaround an issue in OSX: see https://github.com/rstudio/rstudio/issues/6692 .
  cl <- parallel::makeCluster(2, setup_strategy = "sequential")
  parallel::parSapply(cl, 1:2, set.seed)

  unbatched <- simulate(nw0 ~ edges, nsim = 30, coef = -1, control = snctrl(parallel = cl), simplify = FALSE)
  batched1 <- simulate(nw0 ~ edges, nsim = 30, coef = -1, control = snctrl(MCMC.batch = 1, parallel = cl), simplify = FALSE)
  batched6 <- simulate(nw0 ~ edges, nsim = 30, coef = -1, control = snctrl(MCMC.batch = 6, parallel = cl), simplify = FALSE)

  parallel::stopCluster(cl)

  expect_equal(lapply(attr(unbatched, "stats"), dim), lapply(attr(batched1, "stats"), dim))
  expect_equal(lapply(attr(unbatched, "stats"), dim), lapply(attr(batched6, "stats"), dim))

  expect_equal(lengths(unbatched), lengths(batched1))
  expect_equal(lengths(unbatched), lengths(batched6))

  expect_gte(suppressWarnings(ks.test(as.matrix(attr(unbatched, "stats")), as.matrix(attr(batched1, "stats")))$p.value), .05)
  expect_gte(suppressWarnings(ks.test(as.matrix(attr(unbatched, "stats")), as.matrix(attr(batched6, "stats")))$p.value), .05)
})
