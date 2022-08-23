#  File tests/testthat/test-gof.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# Tests of gof()


test_that("gof() works on a bipartite ERGM (#424)", {
  mat <- matrix(c(1,1,0,1,1,0,1,0,0,1,0,1), 4, 3)
  net <-as.network(mat, bipartite=4)
  fit <- ergm(net ~ edges)
  expect_silent(
    r <- gof(fit)
  )
})
