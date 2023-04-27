#  File tests/testthat/test-mple-cov.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

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
  expect_equal(round(coef(summary(m1))[,2], 3), c(0.242, 0.056), ignore_attr = TRUE)
})

test_that("Inverse Hessian from logistic regression model", {
  set.seed(222) # However, this method is not stochastic
  m2 <- ergm(init.sim ~ edges+triangles, estimate = "MPLE",
                control=control.ergm(MPLE.covariance.method = "invHess"))
  expect_equal(round(coef(summary(m2))[,2], 3), c(0.155, 0.034), ignore_attr = TRUE)
})

test_that("Bootstrap covariance method for MPLE", {
  set.seed(333)
  m3 <- ergm(init.sim ~ edges + triangles, estimate = "MPLE",
             control=control.ergm(MPLE.covariance.method = "bootstrap"))
  expect_equal(round(coef(summary(m3))[,2], 3), c(0.257, 0.059), ignore_attr = TRUE)
})
