#  File tests/testthat/test-drop.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

data(sampson)

# Just one covariate. Note that the .mcmc tests mainly test detection
# and overriding of control$force.main. Note that 1/2 has been subtracted from the "maxed" matrices. This is to test detection of non-0-1 extremeness.
test_that("one covariate", {
  samplike.m <- as.matrix(samplike, matrix.type="adjacency")

  maxed.mple <- ergm(samplike~edgecov(samplike.m-1/2))
  expect_equal(coef(maxed.mple), Inf, ignore_attr = TRUE)

  maxed.mcmc <- ergm(samplike~edgecov(samplike.m-1/2),control=control.ergm(force.main=TRUE))
  expect_equal(coef(maxed.mcmc), Inf, ignore_attr = TRUE)

  mined.mple <- ergm(samplike~edgecov(-samplike.m))
  expect_equal(coef(mined.mple), -Inf, ignore_attr = TRUE)

  mined.mcmc <- ergm(samplike~edgecov(-samplike.m),control=control.ergm(force.main=TRUE))
  expect_equal(coef(mined.mcmc), -Inf, ignore_attr = TRUE)
})

# Now, blank out some of the 1s in the matrix so that you still have a
# dropped term, but now multiple parameters are meaningful.
test_that("multiple covariates", {
  samplike.m <- as.matrix(samplike, matrix.type="adjacency")
  samplike.m[4:10,4:10] <- 0

  truth <- c(edges = logit((network.edgecount(samplike)-sum(samplike.m))/(network.dyadcount(samplike)-sum(samplike.m))),
             edgecov.samplike.m = Inf)

  maxed.mple <- ergm(samplike~edges+edgecov(samplike.m))
  expect_equal(coef(maxed.mple), truth)

  maxed.mcmc <- ergm(samplike~edges+edgecov(samplike.m), control=control.ergm(force.main=TRUE, MCMLE.maxit=10))
  expect_equal(coef(maxed.mcmc), truth, tolerance = 0.05)
  expect_equal(logLik(maxed.mcmc), logLik(maxed.mple), tolerance = 0.05, ignore_attr = TRUE)

  truth <- c(edges = logit((network.edgecount(samplike)-sum(samplike.m))/(network.dyadcount(samplike)-sum(samplike.m))),
             `edgecov.-samplike.m` = -Inf)

  mined.mple <- ergm(samplike~edges+edgecov(-samplike.m))
  expect_equal(coef(mined.mple), truth)

  mined.mcmc <- ergm(samplike~edges+edgecov(-samplike.m), control=control.ergm(force.main=TRUE, MCMLE.maxit=10))
  expect_equal(coef(mined.mcmc), truth, tolerance=0.05)
  expect_equal(logLik(mined.mcmc), logLik(mined.mple), tolerance = 0.05, ignore_attr = TRUE)
})

# This is mainly to make sure it doesn't crash for dyad-dependent
# and curved terms.
crash.test <- function() {
  set.seed(1)
  y <- network.initialize(10, directed=FALSE)
  y[1,2]<-y[2,3]<-y[3,4]<-1
  dummy <- ergm(y~edges+triangles+degree(2)+kstar(2)+kstar(5)+gwdegree,
                control=control.ergm(MCMLE.maxit=3), verbose=TRUE) # It doesn't seem to stop for a while.
}

test_that("doesn't crash for dyad-dependent and curved terms", {
  capture_warnings(expect_error(crash.test(), NA))
})
