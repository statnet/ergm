#  File tests/testthat/test-term-project.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
n <- 20
b <- 5

nw0 <- network.initialize(n, bipartite = b, directed = FALSE)
nw0 %v% "b1" <- c(rnorm(b), rep(NA, n-b))
nw0 %v% "b2" <- c(rep(NA, b), rnorm(n-b))
nw1 <- simulate(nw0 ~ edges, coef = 0)
nw2 <- simulate(nw0 ~ edges, coef = 0)
m1 <- as.matrix(nw1)
m2 <- as.matrix(nw2)

nw1.absdiff.b1 <- sum(abs(outer(na.omit(nw1 %v% "b1"), na.omit(nw1 %v% "b1"), FUN = "-")) * tcrossprod(m1))/2
nw1.absdiff.b2 <- sum(abs(outer(na.omit(nw1 %v% "b2"), na.omit(nw1 %v% "b2"), FUN = "-")) * crossprod(m1))/2
nw2.absdiff.b1 <- sum(abs(outer(na.omit(nw2 %v% "b1"), na.omit(nw2 %v% "b1"), FUN = "-")) * tcrossprod(m2))/2
nw2.absdiff.b2 <- sum(abs(outer(na.omit(nw2 %v% "b2"), na.omit(nw2 %v% "b2"), FUN = "-")) * crossprod(m2))/2

f <- nw1 ~ Project(~ absdiff("b1"), 1) + Proj1(~absdiff("b1")) + Project(~absdiff("b2"), 2) + Proj2(~absdiff("b2"))

test_that("Projection(), Proj1(), and Proj2() summary", {
  expect_equal(summary(f),
               c(`Proj1~absdiff.sum.b1` = nw1.absdiff.b1, `Proj1~absdiff.sum.b1` = nw1.absdiff.b1,
                 `Proj2~absdiff.sum.b2` = nw1.absdiff.b2, `Proj2~absdiff.sum.b2` = nw1.absdiff.b2))
})

test_that("Projection(), Proj1(), and Proj2() change", {
  changes <- as.matrix(xor(nw1, nw2), matrix.type = "edgelist")
  expect_equal(c(ergm.godfather(f, changes = changes)),
               rep(c(nw2.absdiff.b1, nw2.absdiff.b2), each = 2))
})
