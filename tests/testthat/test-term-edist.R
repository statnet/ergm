#  File tests/testthat/test-term-edist.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

n <- 50
lat  <- runif(n, 25, 32)
long <- runif(n, 52, 81)

mk_nw <- function(directed = FALSE) {
  nw <- network.initialize(n, directed = directed)
  nw %v% "lat"  <- lat
  nw %v% "long" <- long
  nw
}

test_that("edist summary statistics match Euclidean distance", {
  nw <- san(mk_nw() ~ edges, target.stats = 40)
  el <- as.edgelist(nw)
  d2 <- (lat[el[, 1]] - lat[el[, 2]])^2 + (long[el[, 1]] - long[el[, 2]])^2

  # pow = 1 (default) is the Euclidean distance
  expect_equal(summary(nw ~ edist(c("lat", "long"))), sum(sqrt(d2)), ignore_attr = TRUE)
  # pow = 2 is squared Euclidean distance
  expect_equal(summary(nw ~ edist(c("lat", "long"), pow = 2)), sum(d2), ignore_attr = TRUE)
  # a fractional power
  expect_equal(summary(nw ~ edist(c("lat", "long"), pow = 1.5)), sum(d2^0.75), ignore_attr = TRUE)
  # a vector of powers yields one statistic per power, in order
  expect_equal(summary(nw ~ edist(c("lat", "long"), pow = c(1, 2))),
               c(sum(sqrt(d2)), sum(d2)), ignore_attr = TRUE)
  # coefficient names distinguish the powers
  expect_named(summary(nw ~ edist(c("lat", "long"), pow = c(1, 2))),
               c("edist.lat.long", "edist2.lat.long"))
})

test_that("edist change statistics are correct under add and remove toggles", {
  nw <- mk_nw()
  full <- san(nw ~ edges, target.stats = 40)
  el <- as.edgelist(full)

  # Toggle every edge on, then toggle half of them back off, so the final
  # network is el[-wto, ]; a shuffled order exercises interleaved add/remove.
  wto <- sample(seq_len(NROW(el)), floor(NROW(el) / 2))
  changes <- rbind(el, el[wto, ])
  changes <- changes[sample(seq_len(NROW(changes))), ]

  final <- el[-wto, ]
  d2 <- (lat[final[, 1]] - lat[final[, 2]])^2 + (long[final[, 1]] - long[final[, 2]])^2

  rv <- ergm.godfather(nw ~ edist(c("lat", "long")) + edist(c("lat", "long"), pow = 2),
                       changes = changes)
  expect_equal(as.numeric(rv), c(sum(sqrt(d2)), sum(d2)), ignore_attr = TRUE)
})

test_that("edist requires at least two dimensions", {
  expect_error(summary(mk_nw() ~ edist("lat")), "two or more")
})

test_that("edist supports three or more dimensions", {
  alt <- runif(n, 0, 100)
  nw <- mk_nw()
  nw %v% "alt" <- alt
  full <- san(nw ~ edges, target.stats = 40)
  el <- as.edgelist(full)
  d2 <- (lat[el[, 1]] - lat[el[, 2]])^2 + (long[el[, 1]] - long[el[, 2]])^2 +
    (alt[el[, 1]] - alt[el[, 2]])^2
  expect_equal(summary(full ~ edist(c("lat", "long", "alt"))), sum(sqrt(d2)), ignore_attr = TRUE)
})

test_that("edist works with directed networks", {
  nw <- san(mk_nw(directed = TRUE) ~ edges, target.stats = 40)
  el <- as.edgelist(nw)
  d2 <- (lat[el[, 1]] - lat[el[, 2]])^2 + (long[el[, 1]] - long[el[, 2]])^2
  expect_equal(summary(nw ~ edist(c("lat", "long"))), sum(sqrt(d2)), ignore_attr = TRUE)
})

test_that("edist can be estimated and simulated", {
  nw <- san(mk_nw() ~ edges, target.stats = 60)
  fit <- ergm(nw ~ edges + edist(c("lat", "long")))
  expect_s3_class(fit, "ergm")
  expect_true(all(is.finite(coef(fit))))
  sim <- simulate(fit, nsim = 1)
  expect_gt(network.edgecount(sim), 0)
})
