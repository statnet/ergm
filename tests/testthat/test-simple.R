#  File tests/testthat/test-simple.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# Simulate a network with a high number of nodes with outdegree=3 and a low number with indegree=3:

library(statnet.common)
opttest({

data(sampson)

test_that("extreme outdegree and indegree simulation test", {
  m <- simulate(samplike~odegree(3)+idegree(3), coef=c(100,-100))
  s <- summary(m~odegree(3)+idegree(3))
  expect_lt(diff(s), 0)
})

}, "")
