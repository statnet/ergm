#  File tests/testthat/test-constrain-degrees-edges.R in package ergm, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################


###### Directed
nsim <- 100
n <- 50
m <- 10
d <- 0.1

od <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 1, sum)
id <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 2, sum)
e <- function(nw) network.edgecount(nw)

y0 <- as.network(n, density=d, directed=TRUE)

test_that("is.dyad.independent(ergm) accessor granularity", {
  fit <- suppressWarnings(ergm(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees, estimate="CD",
                               control=control.ergm(CD.maxit=2, CD.nsteps=1, CD.samplesize.per_theta = 10)))
  expect_true(is.dyad.independent(fit, "terms"))
  expect_false(is.dyad.independent(fit, "space"))
  expect_false(is.dyad.independent(fit))
  expect_false(is.na(fit))
  expect_false(anyNA(fit))
})

test_that("constraint = outdegree on a directed network", {
  ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees, coef=setNames(0,""), nsim=nsim, output="stats")
  expect_true(all(sweep(ys[,1:n], 2, od(y0))==0))
  expect_true(any(sweep(ys[,-(1:n)], 2, id(y0))!=0))
})

test_that("constraint = indegree on a directed network", {
  ys <- simulate(y0~receiver(nodes=TRUE)+sender(nodes=TRUE), constraints=~idegrees, coef=setNames(0,""), nsim=nsim, output="stats")
  expect_true(all(sweep(ys[,1:n], 2, id(y0))==0))
  expect_true(any(sweep(ys[,-(1:n)], 2, od(y0))!=0))
})

test_that("both in- and outdegree constraint on a directed network", {
  ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~degrees, coef=setNames(0,""), nsim=nsim, output="stats")
  expect_true(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

  ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees+idegrees, coef=setNames(0,""), nsim=nsim, output="stats")
  expect_true(all(sweep(ys, 2, c(od(y0),id(y0)))==0))
})

test_that("constraint = edges on a directed network", {
  ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~edges, coef=setNames(0,""), nsim=nsim, output="stats")
  # Edges shouldn't vary, but in- and out-degrees should.
  expect_true(all(e(y0)==rowSums(ys[,1:n])))
  expect_true(all(e(y0)==rowSums(ys[,-(1:n)])))
  expect_true(any(sweep(ys[,1:n], 2, od(y0))!=0))
  expect_true(any(sweep(ys[,-(1:n)], 2, id(y0))!=0))
})


test_that("constraint = edges + dyadic constraint on a directed network", {
  y1 <- simulate(y0~edges, coef=0)
  rmat <- matrix(runif(n*n),n,n)
  ys <- simulate(y0~edgecov(!y1) + edgecov(rmat) + sender(nodes=TRUE) + receiver(nodes=TRUE), constraints=~edges + fixallbut(y1), coef=rep(0,n*2+2), nsim=nsim, output="stats")
  # Edges shouldn't vary, but in- and out-degrees should.
  expect_true(all(e(y0)==rowSums(ys[,1:n+2])))
  expect_true(all(e(y0)==rowSums(ys[,1:n+2+n])))
  expect_true(any(sweep(ys[,1:n+2], 2, od(y0))!=0))
  expect_true(any(sweep(ys[,1:n+2+n], 2, id(y0))!=0))
  expect_true(all(ys[,1]==e(y0 & !y1))) # Constrained edges shouldn't vary.
  expect_gt(var(ys[,2]),0) # Unconstrained edges should.
})

###### Undirected
y0 <- as.network(n, density=d, directed=FALSE)

test_that("constraint = degree on an undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~degrees, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(sweep(ys, 2, od(y0))==0))
})

test_that("constraint = edges on an undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~edges, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(e(y0)==rowSums(ys)/2)) # Edges shouldn't vary, but degrees should.
  expect_true(any(sweep(ys, 2, od(y0))!=0))
})

test_that("constraint = edges + dyadic constraint on an undirected network", {
  y1 <- simulate(y0~edges, coef=0)
  rmat <- matrix(runif(n*n),n,n)
  ys <- simulate(y0~edgecov(!y1) + edgecov(rmat) + sociality(nodes=TRUE), constraints=~edges + fixallbut(y1), coef=rep(0,n+2), nsim=nsim, output="stats")
  expect_true(all(e(y0)==rowSums(ys[,-(1:2)])/2)) # Edges shouldn't vary.
  expect_true(any(sweep(ys[,-(1:2)], 2, od(y0))!=0)) # Degrees should vary.
  expect_true(all(ys[,1]==e(y0 & !y1))) # Constrained edges shouldn't vary.
  expect_gt(var(ys[,2]),0) # Unconstrained edges should.
})

###### Bipartite undirected
y0 <- as.network(n, density=d, directed=FALSE, bipartite=m)

test_that("constraint = b1degrees on a bipartite undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b1degrees, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(sweep(ys[,1:m], 2, od(y0))==0))
  expect_true(any(sweep(ys[,-(1:m)], 2, id(y0))!=0))
})

test_that("constraint = b2degrees on a bipartite undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b2degrees, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(sweep(ys[,-(1:m)], 2, id(y0))==0))
  expect_true(any(sweep(ys[,1:m], 2, od(y0))!=0))
})

test_that("both b1 and b2 degrees constraint on a bipartite undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~degrees, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b1degrees+b2degrees, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(sweep(ys, 2, c(od(y0),id(y0)))==0))
})

test_that("constraint = edges on a bipartite undirected network", {
  ys <- simulate(y0~sociality(nodes=TRUE), constraints=~edges, coef=rep(0,n), nsim=nsim, output="stats")
  expect_true(all(e(y0)==rowSums(ys)/2))
  expect_true(any(sweep(ys, 2, c(od(y0),id(y0)))!=0)) # Edges shouldn't vary, but degrees should.
})
