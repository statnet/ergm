#  File tests/testthat/test-term-b12factor.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

bipnet<-network.initialize(4,bipartite=2,directed=FALSE)
add.edges(bipnet,tail=c(1),head=c(4))
# set an attribute on the the vertices of the second mode
set.vertex.attribute(bipnet,'felines',c('cat','tiger'),v=3:4)

test_that("NA check works for b1factor", {
  expect_error(summary(bipnet~b1factor('felines')), "Attribute .\"felines\". has missing data, which is not currently supported by ergm.")
})

test_that("Correct statistics for b2factor()", {
  o <- summary(bipnet~b2factor('felines'))
  expect_equal(o, 1, ignore_attr=TRUE)
})
