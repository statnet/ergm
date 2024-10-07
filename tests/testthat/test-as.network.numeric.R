#  File tests/testthat/test-as.network.numeric.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
test_that("Generating directed network", {
  expect_silent(as.network.numeric(x=5, directed = TRUE))
})
test_that("Generating undirected network", {
  expect_silent(as.network.numeric(x=5, directed = FALSE))
})
test_that("Generating Bipartite network", {
  expect_silent(as.network.numeric(x=5, bipartite = 3, directed = FALSE))
})
test_that("Stop when the loop == TRUE", {
  expect_error(as.network.numeric(x=5, loops = TRUE, bipartite = NULL))
  expect_error(as.network.numeric(x=5, loops = TRUE, directed = TRUE))
  expect_error(as.network.numeric(x=5, loops = TRUE, directed = FALSE))
})
test_that("Stop when the multiple == TRUE", {
  expect_error(as.network.numeric(x=5, multiple = TRUE, bipartite = NULL))
  expect_error(as.network.numeric(x=5, multiple = TRUE, directed = TRUE))
  expect_error(as.network.numeric(x=5, multiple = TRUE, directed = FALSE))
})
test_that("Stop when the hyper == TRUE", {
  expect_error(as.network.numeric(x=5, hyper = TRUE, bipartite = NULL))
  expect_error(as.network.numeric(x=5, hyper = TRUE, directed = TRUE))
  expect_error(as.network.numeric(x=5, hyper = TRUE, directed = FALSE))
})
test_that("Stop when the density is out of [0,1]", {
  expect_error(as.network.numeric(x=5, density = -0.1, bipartite = NULL))
  expect_error(as.network.numeric(x=5, density = -0.1, directed = TRUE))
  expect_error(as.network.numeric(x=5, density = -0.1, directed = FALSE))
  expect_error(as.network.numeric(x=5, density =  1.1, bipartite = NULL))
  expect_error(as.network.numeric(x=5, density =  1.1, directed = TRUE))
  expect_error(as.network.numeric(x=5, density =  1.1, directed = FALSE))
})
test_that("Stop when the number of edges is not integer", {
  expect_error(as.network.numeric(x=5, numedges =  2.5, bipartite = NULL))
  expect_error(as.network.numeric(x=5, numedges =  2.5, directed = TRUE))
  expect_error(as.network.numeric(x=5, numedges =  2.5, directed = FALSE))
})
## # TODO: Add this back when explicit directed && bipartite is an error.
##
## test_that("Stop when a Directed bipartite network is specified by the user", {
##   expect_error(as.network.numeric(x=5, numedges =  6L, directed = TRUE, bipartite = 3))
## })

#test_that("Stop when it is missed by the user to specify whether the network is bipartite or directed", {
#  expect_error(as.network.numeric(x=5, numedges =  6L))
#})

## # TODO: Add this back when warnings can be relied on.
## test_that("Warning message when bipartite == TRUE and user missed to specify whether the network is directed or undirected", {
##   expect_warning(as.network.numeric(x=5, numedges =  6L, bipartite = 3))
## })
test_that("Correct number of edges for bipatite graph",{
  expect_equal(network.edgecount(as.network.numeric(x=5, bipartite = 3, directed = FALSE, numedges = 4)),4)
})
test_that("Correct number of edges for undirected graph",{
  expect_equal(network.edgecount(as.network.numeric(x=5, directed = FALSE, numedges = 10)),10)
})
test_that("Correct number of edges for directed graph",{
  expect_equal(network.edgecount(as.network.numeric(x=5, directed = TRUE, numedges = 18)),18)
})
test_that("Correct number of edges for complete graph or when density == 1", {
  expect_equal(network.edgecount(as.network.numeric(x=5, directed = FALSE, density = 1)),10)
})
test_that("specifying a bigger network that the function can generate results in an error",{
  expect_error(as.network.numeric(x=2^50,directed = FALSE))
})
