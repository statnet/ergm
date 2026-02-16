#  File tests/testthat/test-mple-cov.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

### TODO: Run some very long simulations to get more accurate reference values.

set.seed(14392)
N <- 50
y <- matrix(rbinom(N^2, 1, 0.005), N, N)
diag(y) <- 0
y <- as.network.matrix(y, directed = FALSE)

init.sim <- simulate(y ~ edges + triangles, nsim = 1, coef = c(-log(N) + 4, -0.2))

test_that("Godambe covariance method for MPLE", {
  set.seed(111)
  m1 <- ergm(init.sim ~ edges + triangles, estimate = "MPLE",
              control=control.ergm(MPLE.covariance.method = "Godambe"))
  StdErr1 <- sqrt(diag(vcov(m1)))
  expect_equal(StdErr1, c(0.255, 0.059), ignore_attr = TRUE, tolerance=.05)
})

test_that("Godambe covariance method for MPLE with offset", {
  set.seed(111)
  fit <- ergm(
    init.sim ~ edges + triangles + offset(edges), 
    offset.coef = 0,
    estimate = "MPLE",
    control=control.ergm(MPLE.covariance.method = "Godambe")
  )
  StdErr <- sqrt(diag(vcov(fit)))
  expect_equal(StdErr, c(0.255, 0.059, 0), ignore_attr = TRUE, tolerance=.05)
})

test_that("Inverse Hessian from logistic regression model", {
  set.seed(222) # However, this method is not stochastic
  m2 <- ergm(init.sim ~ edges+triangles, estimate = "MPLE",
                control=control.ergm(MPLE.covariance.method = "invHess"))
  StdErr2 <- sqrt(diag(vcov(m2)))
  expect_equal(StdErr2, c(0.155, 0.034), ignore_attr = TRUE, tolerance=.05)
})

test_that("Bootstrap covariance method for MPLE", {
  set.seed(333)
  m3 <- ergm(init.sim ~ edges + triangles, estimate = "MPLE",
             control=control.ergm(MPLE.covariance.method = "bootstrap"))
  StdErr3 <- sqrt(diag(vcov(m3)))
  expect_equal(StdErr3, c(0.257, 0.060), ignore_attr = TRUE, tolerance=.05)
})

test_that("Bootstrap covariance method for MPLE with offsets", {
  set.seed(445)
  m4 <- ergm(init.sim ~ edges + triangles + offset(triangles), offset.coef=1,
             estimate = "MPLE",
             control=control.ergm(MPLE.covariance.method = "InvHess"))
  StdErr4 <- sqrt(diag(vcov(m4)))
  expect_equal(StdErr4, c(0.155, 0.034, 0), ignore_attr = TRUE, tolerance=.05)
})


