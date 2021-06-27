#  File tests/testthat/test-weighted-population.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

test_that("WtPop produces appropriate samples from uniform random weights", {
  w <- runif(100L)
  n <- 1000000L
  s <- .Call("test_weighted_population", w, n) + 1L
  
  r <- tabulate(s, nbins = length(w))
  
  p <- w/sum(w)
  
  e <- n*p
    
  v <- n*p*(1 - p)
  
  d <- (r - e)/sqrt(v)
  
  expect_true(max(abs(d)) < 6)
})

test_that("WtPop produces appropriate samples from Poisson random weights", {
  w <- as.double(c(rpois(50L, 1), rpois(50L, 10)))
  n <- 1000000L
  s <- .Call("test_weighted_population", w, n) + 1L
  
  r <- tabulate(s, nbins = length(w))
  
  p <- w/sum(w)
  
  e <- n*p
    
  v <- n*p*(1 - p)
  
  v[w == 0] <- 1/100 # so even one sample will fail the test below
  
  d <- (r - e)/sqrt(v)
  
  expect_true(max(abs(d)) < 6)
})
