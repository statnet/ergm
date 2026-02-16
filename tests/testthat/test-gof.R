#  File tests/testthat/test-gof.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

ctrl4 <- control.gof.ergm(nsim=4)

test_that("gof() defaults are correct for bipartite networks", {
  net <- as.network(matrix(c(1,1,0,1,1,0,1,0,0,1,0,1), 4, 3), bipartite=4)
  fit <- ergm(net ~ edges)
  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("b1degree", "b2degree", "dspartners", "distance", "model"))
})



test_that("gof() defaults are correct for undirected networks", {
  data(florentine)
  fit <- ergm((!flomarriage) ~ edges) # Using complement of flomarriage for a dense network test.
  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("degree", "espartners", "distance", "model"))
})



test_that("gof() defaults and GOF handling is correct for directed networks", {
  data(sampson)
  fit <- ergm(samplike ~ edges)

  expect_silent(gof <- gof(fit, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("idegree", "odegree", "espartners", "distance", "model"))

  expect_silent(gof <- gof(fit, GOF=~idegree, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("idegree", "model"))

  expect_silent(gof <- gof(fit, GOF=~idegree-model, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("idegree"))

  expect_silent(gof <- gof(fit, GOF=~model+triadcensus, control=ctrl4))
  expect_silent(plot(gof))
  expect_setequal(which_gof(gof),
                  c("triadcensus", "model"))
})



test_that("gof() is correct for valued networks", {
  data(sampson)
  fit <- ergm(samplike~sum + nonzero + nodematch("group",diff=TRUE,form="sum"),
              response="nominations", reference=~DiscUnif(0, 3),
              control = control.ergm(MCMLE.maxit = 2), eval.loglik = FALSE)

  gof <- gof(fit)
  expect_no_error(expect_no_warning(print(gof)))
  expect_no_error(expect_no_warning(plot(gof)))
  expect_setequal(which_gof(gof), c("model", "cdf"))

  gof <- gof(fit, GOF = ~For(v = -1:4, ~atmost(v)))
  expect_no_error(expect_no_warning(print(gof)))
  expect_no_error(expect_no_warning(plot(gof)))
  expect_setequal(which_gof(gof), c("model", "user"))
})
