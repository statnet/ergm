#  File tests/testthat/test-level-select.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
o <- options(useFancyQuotes=FALSE)

set.seed(123) # Need stable randomization.
data(florentine)
flomarriage %v% "x" <- sample(c(1,11,2,3), 16, replace=TRUE) ## 11 tests for numeric rather than alphabetical sorting.

test_that("Nodal attribute level initialization and sorting", {
  expect_equal(summary(flomarriage ~ nodefactor("x", levels=TRUE)),
               c(nodefactor.x.1 = 3, nodefactor.x.2 = 18, nodefactor.x.3 = 3, nodefactor.x.11 = 16))
})

test_that("Selecting the smallest and largest categories", {
  expect_equal(summary(flomarriage ~ nodefactor("x", levels=SMALLEST)),
               c(nodefactor.x.3 = 3))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=SMALLEST(2))),
               c(nodefactor.x.1 = 3, nodefactor.x.3 = 3))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=LARGEST)),
               c(nodefactor.x.11 = 16))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=LARGEST(2))),
               c(nodefactor.x.2 = 18, nodefactor.x.11 = 16))
})


test_that("Selector negation", {
  expect_equal(summary(flomarriage ~ nodefactor("x", levels=-SMALLEST)),
               c(nodefactor.x.1 = 3, nodefactor.x.2 = 18, nodefactor.x.11 = 16))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=-SMALLEST(2))),
               c(nodefactor.x.2 = 18, nodefactor.x.11 = 16))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=-LARGEST)),
               c(nodefactor.x.1 = 3, nodefactor.x.2 = 18, nodefactor.x.3 = 3))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=-LARGEST(2))),
               c(nodefactor.x.1 = 3, nodefactor.x.3 = 3))
})

test_that("Collapsing categories", {
  expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(2, 5), levels=TRUE)),
               c(nodefactor.x.2 = 18, nodefactor.x.5 = 6, nodefactor.x.11 = 16)) 

  expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(3, 5), levels=TRUE)),
               c(nodefactor.x.5 = 24, nodefactor.x.11 = 16)) 

  expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(2, 5), levels=SMALLEST(2))),
               c(nodefactor.x.2 = 18, nodefactor.x.5 = 6))
})

## Tied categories

set.seed(789) # Need stable randomization.
data(florentine)
flomarriage %v% "x" <- sample(c(1,11,2,3), 16, replace=TRUE) ## 11 tests for numeric rather than alphabetical sorting.

test_that("Tied categories nodal attribute level initialization and sorting", {
  expect_equal(summary(flomarriage ~ nodefactor("x", levels=TRUE)),
               c(nodefactor.x.1 = 2, nodefactor.x.2 = 11, nodefactor.x.3 = 15, nodefactor.x.11 = 12))
})

test_that("Tied categories selecting the smallest and largest", {
  expect_no_warning(expect_equal(summary(flomarriage ~ nodefactor("x", levels=SMALLEST(1))),
                              c(nodefactor.x.1 = 2)))

  expect_warning(expect_equal(summary(flomarriage ~ nodefactor("x", levels=SMALLEST(2))),
                              c(nodefactor.x.1 = 2, nodefactor.x.2 = 11)),
                 "In term 'nodefactor' in package 'ergm': Levels '2' and '11' are tied. Using the order given.")

  expect_no_warning(expect_equal(summary(flomarriage ~ nodefactor("x", levels=SMALLEST(3))),
                                 c(nodefactor.x.1 = 2, nodefactor.x.2 = 11, nodefactor.x.11 = 12)))

  expect_equal(summary(flomarriage ~ nodefactor("x", levels=LARGEST)),
               c(nodefactor.x.3 = 15))

  expect_warning(expect_equal(summary(flomarriage ~ nodefactor("x", levels=LARGEST(2))),
                              c(nodefactor.x.3 = 15, nodefactor.x.11 = 12)),
                 "In term 'nodefactor' in package 'ergm': Levels '2' and '11' are tied. Using the order given.")
})


test_that("Collapsing categories", {
  expect_warning(expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(2, 5), levels=TRUE)),
                              c(nodefactor.x.3 = 15, nodefactor.x.5 = 13, nodefactor.x.11 = 12)),
                 "In term 'nodefactor' in package 'ergm': Levels '2' and '11' are tied. Using the order given.")

  expect_no_warning(expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(3, 5), levels=TRUE)),
                                 c(nodefactor.x.3 = 15, nodefactor.x.5 = 25)))

  expect_warning(expect_warning(expect_equal(summary(flomarriage ~ nodefactor("x" %>% COLLAPSE_SMALLEST(2, 5), levels=SMALLEST(2))),
                                             c(nodefactor.x.3 = 15, nodefactor.x.11 = 12)),
                                "In term 'nodefactor' in package 'ergm': Levels '2' and '11' are tied. Using the order given."),
                 "In term 'nodefactor' in package 'ergm': Levels '3' and '5' are tied. Using the order given.")
})

options(o)
