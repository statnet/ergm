#  File tests/testthat/test-term-diversity.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
diversity <- function(nw, a, which = c("u", "i", "o", "1", "2"), stat = c("range", "distinct")){
  which <- match.arg(which)
  stat <- match.arg(stat)
  stat <- switch(stat,
                 range = function(x) diff(range(x, na.rm=TRUE)),
                 distinct = function(x) length(na.omit(unique(x))))

  m <- as.matrix(nw)
  w <- switch(which,
              u =,
              i =,
              o = nw %v% a,
              `1` = (nw %v% a)[-seq_len(nw%n%"bipartite")],
              `2` = (nw %v% a)[seq_len(nw%n%"bipartite")])

  if(which %in% c("o","1")) m <- t(m)

  r <- suppressWarnings(apply(ifelse(m, w, NA), 2, stat))
  sum(r[!is.infinite(r) & !is.na(r)])
}

data(florentine)

test_that("nodecovrange summary", {
  expect_equal(summary(flomarriage~nodecovrange("wealth")),
               diversity(flomarriage, "wealth", "u", "range"),
               ignore_attr = TRUE)
})

n <- 20
b <- 5

nw0 <- network.initialize(n, bipartite = b, directed = FALSE)
nw0 %v% "b1" <- c(rnorm(b), rep(NA, n-b))
nw0 %v% "b2" <- c(rep(NA, b), rnorm(n-b))
nw1 <- simulate(nw0 ~ edges, coef = 0)

test_that("b1covrange summary", {
  expect_equal(summary(nw1~b1covrange("b2")),
               diversity(nw1, "b2", "1", "range"),
               ignore_attr = TRUE)
})

test_that("b2covrange summary", {
  expect_equal(summary(nw1~b2covrange("b1")),
               diversity(nw1, "b1", "2", "range"),
               ignore_attr = TRUE)
})

data(sampson)
samplike %v% "w" <- rnorm(network.size(samplike))

test_that("nodeocovrange summary", {
  expect_equal(summary(samplike~nodeocovrange("w")),
               diversity(samplike, "w", "o", "range"),
               ignore_attr = TRUE)
})

test_that("nodeicovrange summary", {
  expect_equal(summary(samplike~nodeicovrange("w")),
               diversity(samplike, "w", "i", "range"),
               ignore_attr = TRUE)
})


flomarriage %v% "c" <- sample.int(5, network.size(flomarriage), replace=TRUE)

test_that("nodefactordistinct summary", {
  expect_equal(summary(flomarriage~nodefactordistinct("c")),
               diversity(flomarriage, "c", "u", "distinct"),
               ignore_attr = TRUE)
})

test_that("nodeofactordistinct summary", {
  expect_equal(summary(samplike~nodeofactordistinct("group")),
               diversity(samplike, "group", "o", "distinct"),
               ignore_attr = TRUE)
})

test_that("nodefactordistinct summary", {
  expect_equal(summary(flomarriage~nodefactordistinct("c")),
               diversity(flomarriage, "c", "u", "distinct"),
               ignore_attr = TRUE)
})

test_that("nodeifactordistinct summary", {
  expect_equal(summary(samplike~nodeifactordistinct("group")),
               diversity(samplike, "group", "i", "distinct"),
               ignore_attr = TRUE)
})


n <- 20
b <- 5

nw0 <- network.initialize(n, bipartite = b, directed = FALSE)
nw0 %v% "b1" <- c(sample.int(3, b, TRUE), rep(NA, n-b))
nw0 %v% "b2" <- c(rep(NA, b), sample.int(3, n-b, TRUE))
nw1 <- simulate(nw0 ~ edges, coef = 0)


test_that("b1factordistinct summary", {
  expect_equal(summary(nw1~b1factordistinct("b2")),
               diversity(nw1, "b2", "1", "distinct"),
               ignore_attr = TRUE)
})

test_that("b2factordistinct summary", {
  expect_equal(summary(nw1~b2factordistinct("b1")),
               diversity(nw1, "b1", "2", "distinct"),
               ignore_attr = TRUE)
})
