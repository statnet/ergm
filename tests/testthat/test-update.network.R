#  File tests/testthat/test-update.network.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

aaa <- network.initialize(10,directed=TRUE,loops=TRUE)
aaa %v% 'race' <- rep(c('B','W'),10)
aaa %v% 'letters'<-rep(LETTERS,10)
aaa %v% 'numbers'<-runif(10)
aaa %v% 'astructure'<-rep(list(first="A",second="B"),10)
aaa %v% 'inf.status' <-rbinom(10,1,0.2)
aaa %n% 'verymeta' <- 'so meta'
aaa %n% 'verybeta' <- 'so beta'
aaa[2,1]<-1  # add a single edge

# a random matrix
amat<-structure(c(1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0,  0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1,  0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0,  1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,  0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0), .Dim = c(10L,  10L))

# and edgelist version of the same matrix
ael <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L,  4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L,  7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 10L, 10L,  10L, 10L, 10L, 10L, 1L, 2L, 4L, 5L, 6L, 8L, 3L, 8L, 9L, 6L, 9L,  10L, 1L, 3L, 4L, 6L, 10L, 2L, 3L, 5L, 7L, 8L, 9L, 10L, 3L, 4L,  5L, 7L, 8L, 3L, 4L, 5L, 7L, 8L, 10L, 1L, 4L, 5L, 8L, 9L, 1L,  2L, 1L, 2L, 4L, 5L, 6L, 7L), .Dim = c(48L, 2L))


test_that("update.network() from adjacency matrix", {
  bb<-update(aaa,new = amat,matrix.type = 'adjacency')

  # correct number of edges created
  expect_equal(as.matrix(bb), amat, ignore_attr=TRUE)

  # edges removed
  expect_equal(bb[2,1], 0)

  # original unmodified
  expect_equal(network.edgecount(aaa),1)

  # attributes preserved
  expect_equal(bb%v%'race', aaa%v%'race')
  expect_equal(bb%v%'letters', aaa%v%'letters')
  expect_equal(bb%v%'numbers', aaa%v%'numbers')
  expect_equal(bb%v%'astructure', aaa%v%'astructure')
  expect_equal(bb%n%'verymeta', aaa%n%'verymeta')
  expect_equal(bb%n%'verybeta', aaa%n%'verybeta')

  # flags preserved
  expect_true(has.loops(bb))
})

test_that("update.network() from edgelist", {
  ccc<-update(aaa,ael,matrix.type='edgelist')
  expect_equal(as.matrix(ccc), amat, ignore_attr=TRUE)
})
