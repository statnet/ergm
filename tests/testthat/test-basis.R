#  File tests/testthat/test-basis.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

test_that("basis works as expected", {
  set.seed(0)
  nw <- network(100, directed = FALSE)
  nw %v% "attr" <- rep(1:2, length.out = 100)
  nw2 <- network(100, directed = FALSE)
  nw2 %v% "attr" <- rep(1:2, length.out = 100)
  set.seed(0)
  e <- ergm(nw ~ edges + degree(1) + nodematch("attr"))
  set.seed(0)
  e2 <- ergm(~edges + degree(1) + nodematch("attr"), basis = nw)
  set.seed(0)
  e3 <- ergm(nw2 ~ edges + degree(1) + nodematch("attr"), basis = nw)
  expect_equal(coef(e), coef(e2))
  expect_equal(coef(e2), coef(e3))
  expect_equal(logLik(e), logLik(e2))
  expect_equal(logLik(e2), logLik(e3))
})

test_that("basis works as expected with dyad-dependent constraints", {
  old.opts <- options(ergm.loglik.warn_dyads = FALSE)
  set.seed(0)
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(1:2, length.out = 100)
  nw2 <- network.initialize(100, directed = FALSE)
  nw2 %v% "attr" <- rep(1:2, length.out = 100)
  nw <- san(nw ~ edges, target.stats = c(75), constraints = ~bd(maxout = 3))
  nw2 <- san(nw2 ~ edges, target.stats = c(100), constraints = ~bd(maxout = 3))
  set.seed(0)
  e <- ergm(nw ~ edges + degree(1) + nodematch("attr"), constraints = ~bd(maxout = 3))
  set.seed(0)
  e2 <- ergm(~edges + degree(1) + nodematch("attr"), basis = nw, constraints = ~bd(maxout = 3))
  set.seed(0)
  e3 <- ergm(nw2 ~ edges + degree(1) + nodematch("attr"), basis = nw, constraints = ~bd(maxout = 3))
  expect_equal(coef(e), coef(e2))
  expect_equal(coef(e2), coef(e3))
  expect_equal(logLik(e), logLik(e2))
  expect_equal(logLik(e2), logLik(e3))
  options(old.opts)
})
