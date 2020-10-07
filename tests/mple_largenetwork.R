#  File tests/mple_largenetwork.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
# Note:  n can be made larger if a more stringent test is desired

# First, a test for sparse networks with n edges, so the
# MLE should be log(n / (n*(n-1)-n)) or log(n) - 2*log(n-1)                 
for(n in 500) { # n must be even
  m <- cbind(1:n,n:1) #n must be even here
  e <- ergm(network(m, directed=TRUE)~edges)
  eta <- log(n) - 2*log(n-1)  # What the coefficient estimate should be
  if (round(e$coef - eta, 4) !=0) {
    stop("Failed test of simple MPLE (sparse) estimation for n =",n)
  }
  cat(paste("Passed first MPLE test for n =",n))
}

# Next, a test for less sparse networks in which each
# node has n/100 (out)edges.  Thus the total number of edges is
# n^2/100, so the MLE is log(n^2/100) - log(n(n-1)-n^2/100),
# or 2*log(n) - log(99n^2-100n) = log(n) - log(99n-100)
for(n in 500) { # n should be a multiple of 100
  s <- sample(n-1, n/100)
  m <- cbind(rep(1:n,each=n/100), rep(s, n))
  m[,2] <- m[,2] + apply(m,1,function(x) x[1]<=x[2])
  e <- ergm(network(m, directed=TRUE)~edges)
  eta <- log(n) - log(99*n-100)  # What the coefficient estimate should be
  if (round(e$coef -eta, 4) !=0) {
    stop("Failed test of simple MPLE (nonsparse) estimation for n =",n)
  }
  cat(paste("Passed second MPLE test for n =",n))
}
},"MPLE for large networks")
