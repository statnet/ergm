#  File tests/testthat/test-bridge-target.stats.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

attach(MLE.tools)

library(ergm)
library(statnet.common)

data(florentine)
y <- flomarriage

test_that("Log-likelihood with attainable target statistics",{
  set.seed(1)
  ts <- 3
  llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts)))
  llk <- edges.llk(y,e=ts)
  expect_equal(llk,llk.ergm)
})

test_that("Log-likelihood with unattainable target statistics",{
  set.seed(1)
  ts <- 3.5
  llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts)))
  llk <- edges.llk(y,e=ts)
  expect_equal(llk,llk.ergm,tolerance=0.05)
})

# A nearly empty network with 0 triangles:
nw0 <- network.initialize(10, directed = FALSE)
nw0[1, 2] <- 1

# A network with 8 triangles:
nw1 <- nw0
nw1[cbind(1:9, 2:10)] <- 1
nw1[cbind(1:8, 3:10)] <- 1

# mle <-coef( ergm(nw1~triangles, eval.loglik = FALSE))
mle <- -0.2144383 # hard-code to save time

set.seed(1)

for (theta in c(mle, rnorm(1, -0.5, 0.25))) { # MLE and a random value
  test_that(paste("log-likelihood calculation for the situation where network stats differ significantly from target stats:",
                  if (theta == mle) "MLE" else format(theta)), {
    # network stats != target stats
    llk0 <- ergm.bridge.dindstart.llk(nw0~triangles, coef = theta, target.stats = 8, llkonly = FALSE)
    # network stats == target stats
    llk1 <- ergm.bridge.dindstart.llk(nw1~triangles, coef = theta, llkonly = FALSE)

    # difference |Z| < 3
    expect_lt((llk0$llk - llk1$llk)^2 / (llk0$vcov.llr + llk1$vcov.llr), 9)
  })
}

detach(MLE.tools)
