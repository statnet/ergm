#  File tests/testthat/test-term-undirected.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

# an undirected nw
data(faux.mesa.high)
fmh <- faux.mesa.high
set.seed(7)
set.edge.attribute(fmh, "GradeMet", rbinom(203, 6, .5))

# a small undirected nw w/ lots o' triangles
set.seed(20)
t<-trunc(runif(160, 1, 20))
set.seed(21)
h<-trunc(runif(160, 1, 20))
el <- cbind(t,h)
bad <- which(el[,2]==el[,1])
el[bad,2] = el[bad,2]+1
unnw <- network(el, directed=FALSE)
unnw %v% "Pet" <- c("dog", "cat")

test_that("altkstar, undirected, ", {
  s.0 <- gf_summary(fmh~altkstar)
  e.0 <- ergm(fmh~altkstar(1, fixed=TRUE), estimate="MPLE")
  e.l <- ergm(fmh~altkstar(.5, fixed=TRUE), estimate="MPLE")
  s.f <- gf_summary(fmh~altkstar(1, fixed=TRUE))
  e.lf <- ergm(fmh~altkstar(.9, fixed=TRUE), estimate="MPLE")

  expect_summary(s.0[1:10], e.0, c(51,30,28,18,10,2,4,1,2,1), -3.234)
  expect_equal(coef(e.l), -4.166, tolerance=0.001, ignore_attr=TRUE)
  expect_equal(s.f, 258, ignore_attr=TRUE)
  expect_equal(coef(e.lf), -3.494, tolerance=0.001, ignore_attr=TRUE)
})

test_that("concurrent, undirected", {
  s.0 <- gf_summary(fmh~concurrent)
  e.0 <- ergm(fmh~concurrent, estimate="MPLE")
  s.b <- gf_summary(fmh~concurrent(by=function(x) x %v% "Grade"))
  e.b <- ergm(fmh~concurrent(by="Sex"), estimate="MPLE")

  expect_summary(s.0, e.0, 97, -4.871)
  expect_summary(s.b, e.b, c(35,15,18,8,13,8), -c(5.17301, 4.67697))
})

test_that("concurrentties, undirected", {
  s.0 <- gf_summary(fmh~concurrentties)
  e.0 <- ergm(fmh~concurrentties, estimate="MPLE")
  s.b <- gf_summary(fmh~concurrentties(by="Grade"))
  e.b <- ergm(fmh~concurrentties(by=~Sex), estimate="MPLE")

  expect_summary(s.0, e.0, 258, -3.234)
  expect_summary(s.b, e.b, c(103,51,36,19,31,18), -c(3.078,3.429))
})

test_that("cyclicalties, directed", {
  s.0 <- gf_summary(fmh~cyclicalties)
  e.0 <- ergm(fmh~cyclicalties, estimate="MPLE")
  s.a <- gf_summary(fmh~cyclicalties("Race"))
  e.a <- ergm(fmh~cyclicalties("Race"), estimate="MPLE")

  expect_summary(s.0, e.0, 120, -0.4868)
  expect_summary(s.a, e.a, 40, -0.4430)
})

test_that("degree, undirected", {
  s.d <- gf_summary(fmh~degree(2:3))
  e.d <- ergm(fmh~degree(0), estimate="MPLE")
  s.db <- gf_summary(fmh~degree(1:3, function(x) x %v% "Grade"))
  e.db <- ergm(fmh~degree(4, "Sex"), estimate="MPLE")
  s.dbh <- gf_summary(fmh~degree(4:5, by="Sex", homophily=TRUE))
  e.dbh <- ergm(fmh~degree(2, by=~Grade, homophily=TRUE), estimate="MPLE")

  expect_summary(s.d, e.d, c(30,28), 5.11)
  expect_summary(s.db, e.db, c(15,9,9,9,4,2,11,5,9,9,4,2,5,5,4,2,3,2), -c(.345, .6005))
  expect_summary(s.dbh, e.dbh, c(10,3), -.5713)
})

test_that("degrange, undirected", {
  s.0 <- gf_summary(fmh~degrange(1:3))
  e.0 <- ergm(fmh~degrange(1:3), estimate="MPLE")
  s.h <- gf_summary(fmh~degrange(1:3, by="Sex", homophily=TRUE))
  e.h <- ergm(fmh~degrange(1:3, by=~Sex, homophily=TRUE), estimate="MPLE")

  expect_summary(s.0, e.0, c(148, 97, 67), -c(4.349, 4.067, 3.178  ))
  expect_summary(s.h, e.h, c(122, 65, 36), -c(3.389, 3.032, 2.368 ))
})

test_that("degcrossprod, undirected", {
  s.0 <- gf_summary(unnw~degcrossprod)
  e.0 <- ergm(unnw~degcrossprod, estimate="MPLE")

  expect_summary(s.0, e.0, c(56.30102), c(0.099))
})

test_that("degcor, undirected", {
  s.0 <- gf_summary(unnw~degcor)
  e.0 <- ergm(unnw~degcor, estimate="MPLE")

  expect_summary(s.0, e.0, -c(0.09789041 ), c(2.337))
})

test_that("degcor handles an empty initial network", {
  empty <- network.initialize(5, directed=FALSE)
  mple <- ergmMPLE(empty ~ edges + degcor, output="matrix")

  expect_true(all(is.finite(mple$predictor)))
  expect_equal(as.numeric(summary(empty ~ degcor)), 0)
})

test_that("degcor change statistics account for variance changes", {
  nw <- network.initialize(6, directed=FALSE)
  network::add.edges(nw, c(1, 1, 1, 2), c(2, 3, 4, 3))
  mple <- ergmMPLE(nw ~ edges + degcor, output="array")
  statistic <- as.numeric(summary(nw ~ degcor))

  added <- network::network.copy(nw)
  network::add.edges(added, 1, 5)
  expect_equal(mple$predictor[1, 5, "degcor"],
               as.numeric(summary(added ~ degcor)) - statistic)

  removed <- network::network.copy(nw)
  network::delete.edges(removed, network::get.edgeIDs(removed, v=1, alter=2))
  expect_equal(mple$predictor[1, 2, "degcor"],
               statistic - as.numeric(summary(removed ~ degcor)))
})

test_that("degcor updates its variance during sampling", {
  nw <- network.initialize(6, directed=FALSE)
  network::add.edges(nw, c(1, 1, 1, 2), c(2, 3, 4, 3))
  set.seed(1)
  sim <- simulate(nw ~ edges + degcor, coef=c(0, 0), nsim=5, output="network",
                  control=control.simulate.formula(MCMC.burnin=20, MCMC.interval=2))

  expect_true(all(vapply(sim, function(y) is.finite(as.numeric(summary(y ~ degcor))), logical(1))))
})

test_that("degree1.5, undirected", {
  s.0 <- gf_summary(fmh~degree1.5)
  e.0 <- ergm(fmh~degree1.5, estimate="MPLE")

  expect_summary(s.0, e.0, 795.7458, -1.1398)
})

test_that("gwdegree, undirected", {
  s.d <- gf_summary(fmh~gwdegree())
  expect_error(summary(fmh~gwdegree(cutoff=9)), ".*Term .gwdegree. has encountered a network for which degree of some node exceeded the cut-off setting of 9. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  e.d <- ergm(fmh~gwdegree(.4, fixed=TRUE), estimate="MPLE")
  s.df <- gf_summary(fmh~gwdegree(.3, fixed=TRUE))
  e.df <- ergm(fmh~gwdegree(.2, fixed=TRUE), estimate="MPLE")
  s.dfa <- gf_summary(fmh~gwdegree(.1, fixed=TRUE, attr=function(x) x %v% "Grade"))
  e.dfa <- ergm(fmh~gwdegree(.1, fixed=TRUE, attr=~Grade), estimate="MPLE")

  expect_summary(head(s.d), e.d, setNames(c(51,30,28,18,10,2), paste0("gwdegree#",1:6)), c(gwdeg.fixed.0.4=-13.59067))
  expect_summary(s.df, e.df, 178.4312, -18.2508)
  expect_summary(s.dfa, e.dfa,
    c(53.58148, 25.53534, 30.83418, 17.79934, 19.31326, 10.80933),
    -c(23.94060, 23.30646, 23.51430, 23.31140, 25.11103, 26.88088))

  # Also, check that gwdegree(0) is evaluates correctly.
  expect_equal(summary(fmh ~ gwdegree(0, fixed=TRUE)), network.size(fmh) - summary(fmh ~ degree(0)), ignore_attr = TRUE)
})

test_that("kstar, undirected", {
  s.k <- gf_summary(fmh~kstar(1:3))
  e.k <- ergm(fmh~kstar(c(2,4)), estimate="MPLE")
  s.ka <- gf_summary(fmh~kstar(2, "Grade"))
  e.ka <- ergm(fmh~kstar(2, "Sex"), estimate="MPLE")

  expect_summary(s.k, e.k, c(406, 659, 1010), c(-1.45086, .06255))
  expect_summary(s.ka, e.ka, 466, -1.535175)
})

test_that("opentriad, undirected", {
  s.0 <- gf_summary(fmh~opentriad)
  e.0 <- ergm(fmh~opentriad, estimate="MPLE")

  expect_summary(s.0, e.0, 473, -1.652)
})

test_that("sociality, undirected", {
  s.0 <- gf_summary(fmh~sociality)
  s.b <- gf_summary(fmh~sociality(nodes=-(2:203)))

  expect_equal(head(s.0), c(4,0,0,1,0,0), ignore_attr=TRUE)
  expect_equal(s.b, c(13,3,1), ignore_attr=TRUE)
})

test_that("transitiveties, directed", {
  s.0 <- gf_summary(fmh~transitiveties)
  e.0 <- ergm(fmh~transitiveties, estimate="MPLE")
  s.a <- gf_summary(fmh~transitiveties("Race"))
  e.a <- ergm(fmh~transitiveties("Race"), estimate="MPLE")

  expect_summary(s.0, e.0, 120, -0.4868)
  expect_summary(s.a, e.a, 40, -0.4430)
})

test_that("tripercent, undirected", {
  s.0 <- gf_summary(unnw~tripercent)
  e.0 <- ergm(unnw~tripercent, estimate="MPLE")
  s.a <- gf_summary(unnw~tripercent("Pet"))
  e.a <- ergm(unnw~tripercent(~Pet), estimate="MPLE")

  expect_summary(s.0, e.0, 29.19463, 0.4492)
  expect_summary(s.a, e.a, 29.09091, 0.2501)
})
