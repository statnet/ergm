#  File tests/testthat/test-ergm.bridge.llr.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

attach(MLE.tools)

n <- 10
nw <- network.initialize(n)
nw[cbind(1:(n-1), 2:n)] <- 1

test_that("log-likelihood ratio for an edges model from 0 to MLE", {
  expect_equal(
    ergm.bridge.llr(nw ~ edges, from = 0, to = edges.theta(nw))$llr,
    edges.llk(nw) - edges.llk(nw, 0),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

test_that("log-likelihood ratio for an edges model between two random values", {
  expect_equal(
    ergm.bridge.llr(nw ~ edges, from = (from <- rnorm(1)), to = (to <- rnorm(1)))$llr,
    edges.llk(nw, to) - edges.llk(nw, from),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

test_that("log-likelihood for an edges model MLE with missing data", {
  expect_equal(
    ergm.bridge.dindstart.llk(nw~edges, coef=edges.theta(nw)),
    edges.llk(nw),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

nw[cbind(2:n, 1:(n-1))] <- NA

test_that("log-likelihood ratio for an edges model from 0 to MLE with missing data", {
  expect_equal(
    ergm.bridge.llr(nw ~ edges, from = 0, to = edges.theta(nw))$llr,
    edges.llk(nw) - edges.llk(nw, 0),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

test_that("log-likelihood ratio for an edges model between two random values with missing data", {
  expect_equal(
    ergm.bridge.llr(nw ~ edges, from = (from <- rnorm(1)), to = (to <- rnorm(1)))$llr,
    edges.llk(nw, to) - edges.llk(nw, from),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

test_that("log-likelihood for an edges model MLE with missing data", {
  expect_equal(
    ergm.bridge.dindstart.llk(nw~edges, coef=edges.theta(nw)),
    edges.llk(nw),
    ignore_attr = TRUE, tolerance = 0.01
  )
})

detach(MLE.tools)
