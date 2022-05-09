#  File tests/testthat/test-term-edgecov.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

test_that("edgecov works with undirected unipartite networks", {
  nw <- network.initialize(5, dir=FALSE)
  changes <- matrix(c(1,2,
                      2,3,
                      4,5,
                      1,2), ncol = 2, byrow = TRUE)
  m <- matrix(0L, 5, 5)
  m[upper.tri(m)] <- 1:10
  m <- m + t(m) # documentation says undirected case assumes covariate is undirected
  rv <- ergm.godfather(nw ~ edgecov(m), changes=changes)
  expect_identical(as.integer(rv), sum(m[changes[c(2,3),]]))
})

test_that("edgecov works with undirected bipartite networks", {
  nw <- network.initialize(10, bip=4, dir=FALSE)
  changes <- matrix(c(1,5,
                      1,9,
                      3,5,
                      2,6,
                      4,10,
                      4,8,
                      3,7,
                      1,9,
                      3,5,
                      3,6), ncol = 2, byrow = TRUE)
  m <- matrix(c(-(1:12), 1:12), nrow=4, ncol=6)
  rv <- ergm.godfather(nw ~ edgecov(m), changes=changes)
  rtk <- c(1,4,5,6,7,10)
  expect_identical(as.integer(rv), sum(m[cbind(changes[rtk,1], changes[rtk,2] - 4)]))
})

test_that("edgecov works with directed networks", {
  nw <- network.initialize(5, dir=TRUE)
  changes <- matrix(c(1,2,
                      3,2,
                      4,5,
                      5,1,
                      2,3,
                      1,2,
                      5,1), ncol = 2, byrow = TRUE)
  m <- matrix(0L, 5, 5)
  m[upper.tri(m)] <- 1:10
  m[lower.tri(m)] <- -(1:10)
  rv <- ergm.godfather(nw ~ edgecov(m), changes=changes)
  expect_identical(as.integer(rv), sum(m[changes[c(2,3,5),]]))
})
