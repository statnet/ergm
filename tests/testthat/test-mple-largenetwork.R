#  File tests/testthat/test-mple-largenetwork.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

# Note:  n can be made larger if a more stringent test is desired

# First, a test for sparse networks with n edges, so the
# MLE should be log(n / (n*(n-1)-n)) or log(n) - 2*log(n-1)
for(n in 500) { # n must be even
  test_that(paste0("first MPLE test for n =", n), {
    m <- cbind(1:n,n:1) #n must be even here
    e <- ergm(network(m, directed=TRUE)~edges)
    eta <- log(n) - 2*log(n-1)  # What the coefficient estimate should be
    expect_equal(coef(e), eta, tolerance=0.0001, ignore_attr=TRUE)
  })
}

# Next, a test for less sparse networks in which each
# node has n/100 (out)edges.  Thus the total number of edges is
# n^2/100, so the MLE is log(n^2/100) - log(n(n-1)-n^2/100),
# or 2*log(n) - log(99n^2-100n) = log(n) - log(99n-100)
for(n in 500) { # n should be a multiple of 100
  test_that(paste0("second MPLE test for n =", n), {
    s <- sample(n-1, n/100)
    m <- cbind(rep(1:n,each=n/100), rep(s, n))
    m[,2] <- m[,2] + apply(m,1,function(x) x[1]<=x[2])
    e <- ergm(network(m, directed=TRUE)~edges)
    eta <- log(n) - log(99*n-100)  # What the coefficient estimate should be
    expect_equal(coef(e), eta, tolerance=0.0001, ignore_attr=TRUE)
  })
}
