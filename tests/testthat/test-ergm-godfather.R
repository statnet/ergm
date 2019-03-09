#  File tests/testthat/test-ergm-godfather.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
context("test-ergm-godfather.R")

nw <- network.initialize(4, dir=TRUE)
nw[,,names.eval="w",add.edges=TRUE] <-
  matrix(c(0,1,2,0,
           0,0,1,3,
           1,0,0,0,
           1,4,1,0),
         4,4,
         byrow=TRUE)

test_that("ergm.godfather() with toggles", {
  gf <- ergm.godfather(nw~edges+triangles, changes=list(matrix(c(1,2,
                                                                 2,3),
                                                               ncol=2,byrow=TRUE),
                                                        matrix(c(1,3,
                                                                 2,3),
                                                               ncol=2,byrow=TRUE)),
                       stats.start=TRUE)

  expect_equivalent(as.matrix(gf),
                    matrix(c(8,8,
                             6,2,
                             6,3),
                           ncol=2,byrow=TRUE))
})


test_that("ergm.godfather() with values", {
  gf <- ergm.godfather(nw~edges+triangles, changes=list(matrix(c(1,2,1,
                                                                 2,3,0),
                                                               ncol=3,byrow=TRUE),
                                                        matrix(c(1,3,1,
                                                                 2,3,1),
                                                               ncol=3,byrow=TRUE)))

  expect_equivalent(as.matrix(gf),
                    matrix(c(7,4,
                             8,8),
                           ncol=2,byrow=TRUE))
})


test_that("valued ergm.godfather()", {
  gf <- ergm.godfather(nw~nonzero+sum+transitiveweights(), response="w",
                       changes=list(matrix(c(1,2,1,
                                             2,3,0),
                                           ncol=3,byrow=TRUE),
                                    matrix(c(1,3,1,
                                             2,3,1),
                                           ncol=3,byrow=TRUE)))

  expect_equivalent(as.matrix(gf),
                    matrix(c(7,13,3,
                             8,13,5),
                           ncol=3,byrow=TRUE))
})


test_that("valued ergm.godfather() returning the network", {
  gf <- ergm.godfather(nw~nonzero+sum+transitiveweights(), response="w",
                       changes=list(matrix(c(1,2,1,
                                             2,3,0),
                                           ncol=3,byrow=TRUE),
                                    matrix(c(1,3,1,
                                             2,3,1),
                                           ncol=3,byrow=TRUE)),
                       end.network=TRUE)

  expect_equivalent(as.matrix(gf, attrname="w"),
                    matrix(c(0,1,1,0,
                             0,0,1,3,
                             1,0,0,0,
                             1,4,1,0),
                           4,4,
                           byrow=TRUE)
                    )
})

