#  File tests/testthat/test-gof.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
# Tests of gof()

ctrl4 <- control.gof.ergm(nsim=4)

test_that("gof() defaults are correct for bipartite networks", {
  net <- as.network(matrix(c(1,1,0,1,1,0,1,0,0,1,0,1), 4, 3), bipartite=4)
  fit <- ergm(net ~ edges)
  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(as.character(statnet.common::list_rhs.formula(gof$GOF)),
                  c("b1degree", "b2degree", "espartners", "distance", "model"))
})



test_that("gof() defaults are correct for undirected networks", {
  data(florentine)
  fit <- ergm((!flomarriage) ~ edges) # Using complement of flomarriage for a dense network test.
  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(as.character(statnet.common::list_rhs.formula(gof$GOF)),
                  c("degree", "espartners", "distance", "model"))
})



test_that("gof() defaults and GOF handling is correct for directed networks", {
  data(sampson)
  fit <- ergm(samplike ~ edges)

  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(as.character(list_rhs.formula(gof$GOF)),
                  c("idegree", "odegree", "espartners", "distance", "model"))

  expect_silent(gof <- gof(fit, GOF=~idegree, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(as.character(list_rhs.formula(gof$GOF)),
                  c("idegree", "model"))

  expect_silent(gof <- gof(fit, GOF=~idegree-model, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(as.character(list_rhs.formula(gof$GOF)),
                  c("idegree"))
})
