#  File tests/testthat/test-term-directed.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

# a directed nw
load("sampson.wrong.RData") # Old (wrong) version of sampson's monks
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")

test_that("asymmetric, directed", {
  s.0 <- summary(samplike~asymmetric)
  e.0 <- ergm(samplike~asymmetric, estimate="MPLE")
  s.a <- summary(samplike~asymmetric("group"))
  e.a <- ergm(samplike~asymmetric("group"), estimate="MPLE")
  s.ad <- summary(samplike~asymmetric(function(x) x %v% "group", diff=TRUE))
  e.ad <- ergm(samplike~asymmetric(~group, diff=TRUE), estimate="MPLE")
  s.ak <- summary(samplike~asymmetric("group", levels=3))
  e.ak <- ergm(samplike~asymmetric(function(x) x %v% "group", levels=3), estimate="MPLE")
  s.adk <- summary(samplike~asymmetric(~group, diff=TRUE, keep=1:2))
  e.adk <- ergm(samplike~asymmetric("group", diff=TRUE, keep=c(1,3)), estimate="MPLE")
  expect_summary(s.0, e.0, 32, -1.33)
  expect_summary(s.a, e.a, 17, -0.6008)
  expect_summary(s.ad, e.ad, c(7,2,8), c(-.6931, -.6931, -.4855))
  expect_summary(s.ak, e.ak, 8, -.4855)
  expect_summary(s.adk, e.adk, c(7,2), c(-.6931, -.4855))
})

test_that("ctripe=ctriad, directed", {
  s.0 <- summary(samplike~ctriple)
  e.0 <- ergm(samplike~ctriad, estimate="MPLE")
  s.a <- summary(samplike~ctriple(function(x) x %v% "group"))
  e.a <- ergm(samplike~ctriple("group"), estimate="MPLE")
  s.ad <- summary(samplike~ctriad(~group, diff=TRUE))
  e.ad <- ergm(samplike~ctriple(function(x) x %v% "group", diff=TRUE), estimate="MPLE")
  expect_summary(s.0, e.0, 39, -.3522)
  expect_summary(s.a, e.a, 34, .1217)
  expect_summary(s.ad, e.ad, c(8,4,22), -c(.1949, -.6931, -.2023))
})

test_that("cyclicalties, directed", {
  s.0 <- summary(samplike~cyclicalties)
  e.0 <- ergm(samplike~cyclicalties, estimate="MPLE")
  s.b <- summary(samplike~cyclicalties(attr=function(x) x %v% "group"))
  e.b <- ergm(samplike~cyclicalties(attr="group"), estimate="MPLE")
  expect_summary(s.0, e.0, 62, -.4154)
  expect_summary(s.b, e.b, 55, .2289)
})

test_that("idegrange, directed", {
  s.0 <- summary(samplike~idegrange(5:8))
  e.0 <- ergm(samplike~idegrange(5:8), estimate="MPLE")
  s.h <- summary(samplike~idegrange(5:8, by="group", homophily=TRUE))
  e.h <- ergm(samplike~idegrange(5:8, by=~group, homophily=TRUE), estimate="MPLE")
  expect_summary(s.0, e.0, c(9, 6, 4, 3), -c(-0.1431, 1.0986, 1.1451, 0.2231 ))
  expect_summary(s.h, e.h, c(5, 3, 0, 0), -c(-0.5108, -2.1972 , Inf, Inf ))
})

test_that("gwidegree, directed", {
  expect_silent(s.d <- summary(samplike~gwidegree(cutoff=11)))
  expect_error(summary(samplike~gwidegree(cutoff=7)), ".*Term .gwidegree. has encountered a network for which in-degree of some node exceeded the cut-off setting of 7. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  e.d <- ergm(samplike~gwidegree(.4, fixed=TRUE), estimate="MPLE")
  s.df <- summary(samplike~gwidegree(.3, fixed=TRUE))
  e.df <- ergm(samplike~gwidegree(.2, fixed=TRUE), estimate="MPLE")
  s.dfa <- summary(samplike~gwidegree(.1, TRUE, "group"))
  e.dfa <- ergm(samplike~gwidegree(.5, TRUE, function(x) x %v% "group"), estimate="MPLE")
  expect_summary(head(s.d), e.d, setNames(c(0,3,5,1,3,2), paste0("gwidegree#",1:6)), c(gwideg.fixed.0.4=-5.783202))
  expect_summary(s.df, e.df, 23.89614, -2.247936)
  expect_summary(s.dfa, e.dfa, c(7.715119, 4.408762, 7.734290), -c(5.460448, 5.754111, 6.144961))
})

test_that("gwodegree, directed", {
  s.d <- summary(samplike~gwodegree())
  expect_error(summary(samplike~gwodegree(cutoff=5)), ".*Term .gwodegree. has encountered a network for which out-degree of some node exceeded the cut-off setting of 5. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  e.d <- ergm(samplike~gwodegree(.4, fixed=TRUE), estimate="MPLE")
  s.df <- summary(samplike~gwodegree(.3, fixed=TRUE))
  e.df <- ergm(samplike~gwodegree(.2, fixed=TRUE), estimate="MPLE")
  s.dfa <- summary(samplike~gwodegree(.1, TRUE, ~group))
  e.dfa <- ergm(samplike~gwodegree(.5, TRUE, "group"), estimate="MPLE")
  expect_summary(head(s.d), e.d, setNames(c(0,0,1,5,7,5), paste0("gwodegree#",1:6)), c(gwodeg.fixed.0.4=-1.990492))
  expect_summary(s.df, e.df, 24.23040, 43.61801)
  expect_summary(s.dfa, e.dfa, c(7.735906, 4.419631, 7.736070), -c(4.1860720, 5.9706455, 0.4921623))

  # Also, check that gwodegree(0) works correctly, including when there are none.
  expect_equal(summary(samplike~gwodegree(0, fixed = TRUE)), network.size(samplike), ignore_attr = TRUE)
  expect_equal(coef(ergm(samplike~edges + gwodegree(0, fixed = TRUE), estimate="MPLE"))[[2]], +Inf)
})

test_that("idegree, directed", {
  s.d <- summary(samplike~idegree(2:3))
  e.d <- ergm(samplike~idegree(2), estimate="MPLE")
  s.db <- summary(samplike~idegree(1:3, "group"))
  (e.db <- ergm(samplike~idegree(3, function(x) x %v% "group"), estimate="MPLE")) |>
    expect_warning("The MPLE does not exist!")
  s.dbh <- summary(samplike~idegree(4:5, "group", TRUE))
  e.dbh <- ergm(samplike~idegree(2, ~group, TRUE), estimate="MPLE")
  expect_summary(s.d, e.d, c(3,5), 1.223775)
  expect_summary(s.db, e.db, c(0,2,1,0,1,2,0,0,2), c(-0.6931472, 0.8183103, 17.4324836))
  expect_summary(s.dbh, e.dbh, c(3,2), .3448)
})

test_that("intransitive, directed", {
  s.0 <- summary(samplike~intransitive)
  e.0 <- ergm(samplike~intransitive, estimate="MPLE")
  expect_summary(s.0, e.0, 224, -.21)
})

test_that("idegree1.5, directed", {
  s.0 <- summary(samplike~idegree1.5)
  e.0 <- ergm(samplike~idegree1.5, estimate="MPLE")
  expect_summary(s.0, e.0, 214.6543, -.2387)
})

test_that("istar, directed", {
  s.k <- summary(samplike~istar(1:3))
  e.k <- ergm(samplike~istar(c(2,4)), estimate="MPLE")
  s.ka <- summary(samplike~istar(2, function(x) x %v% "group"))
  e.ka <- ergm(samplike~istar(2, "group"), estimate="MPLE")
  expect_summary(s.k, e.k, c(88,233,455), c(-.28615, .02477))
  expect_summary(s.ka, e.ka, 100, .2401)
})

test_that("m2star, directed", {
  s.0 <- summary(samplike~m2star)
  e.0 <- ergm(samplike~m2star, estimate="MPLE")
  expect_summary(s.0, e.0, 378, -.1028)
})

test_that("mutual, directed", {
  s.0 <- summary(samplike~mutual)
  e.0 <- ergm(samplike~mutual, estimate="MPLE")
  s.s <- summary(samplike~mutual(same=function(x) x %v% "group"))
  e.s <- ergm(samplike~mutual(same="group"), estimate="MPLE")
  s.b <- summary(samplike~mutual(by=~Trinity))
  e.b <- ergm(samplike~mutual(by="Trinity"), estimate="MPLE")
  s.sd <- summary(samplike~mutual(same="group", diff=TRUE))
  e.sd <- ergm(samplike~mutual(same=function(x) x %v% "group", diff=TRUE), estimate="MPLE")
  s.sk <- summary(samplike~mutual(same="group", levels=2))
  e.sk <- ergm(samplike~mutual(same=~group, levels=1), estimate="MPLE")
  s.bk <- summary(samplike~mutual(by="Trinity", keep=2))
  e.bk <- ergm(samplike~mutual(by="Trinity", keep=2:3), estimate="MPLE")
  expect_summary(s.0, e.0, 28, .5596)
  expect_summary(s.s, e.s, 23, .9954)
  expect_summary(s.b, e.b, c(17,18,21), c(.0157, .0667, .8077))
  expect_summary(s.sd, e.sd, c(8,4,11), c(.8267, 1.3863, 1.0116))
  expect_summary(s.sk, e.sk, 4, .8266)
  expect_summary(s.bk, e.bk, 18, c(.0714, .8108))
})

test_that("nearsimmelian, directed", {
  s.0 <- summary(samplike~nearsimmelian)
  e.0 <- ergm(samplike~nearsimmelian, estimate="MPLE")
  expect_summary(s.0, e.0, 18, -.4366483)
})

test_that("nodeicov, directed", {
  s.a <- summary(samplike~nodeicov("YearsServed"))
  e.a <- ergm(samplike~nodeicov(function(x) x %v% "YearsServed"), estimate="MPLE")
  s.at <- summary(samplike~nodeicov(~YearsServed^2))
  e.at <- ergm(samplike~nodeicov(~(.%v%"YearsServed")^2), estimate="MPLE")
  s.att <- summary(samplike~nodeicov(~poly(YearsServed,2,raw=TRUE)))
  expect_summary(s.a, e.a, 439, -.1739)
  expect_summary(s.at, e.at, 2345, -.02805)
  expect_equal(s.att, c(439,2345), ignore_attr=TRUE)
})

test_that("nodeifactor, directed", {
  s.a <- summary(samplike~nodeifactor("group"))
  e.a <- ergm(samplike~nodeifactor(~group), estimate="MPLE")
  s.ab <- summary(samplike~nodeifactor(function(x) x %v% "Trinity", levels=TRUE))
  e.ab <- ergm(samplike~nodeifactor("Trinity", base=(2:3)), estimate="MPLE")
  expect_summary(s.a, e.a, c(13, 46), -c(1.4424, .4618))
  expect_summary(s.ab, e.ab, c(28, 29, 31), -.9719)
})

test_that("nodeocov, directed", {
  s.a <- summary(samplike~nodeocov("YearsServed"))
  e.a <- ergm(samplike~nodeocov(function(x) x %v% "YearsServed"), estimate="MPLE")
  s.at <- summary(samplike~nodeocov(~YearsServed^2))
  e.at <- ergm(samplike~nodeocov(~(.%v%"YearsServed")^2), estimate="MPLE")
  s.att <- summary(samplike~nodeocov(~poly(YearsServed,2,raw=TRUE)))
  expect_summary(s.a, e.a, 467, -.1581)
  expect_summary(s.at, e.at, 2691, -.02243)
  expect_equal(s.att, c(467, 2691), ignore_attr=TRUE)
})

test_that("nodeofactor, directed", {
  s.a <- summary(samplike~nodeofactor("group"))
  e.a <- ergm(samplike~nodeofactor(~group), estimate="MPLE")
  s.ab <- summary(samplike~nodeofactor("Trinity", levels=TRUE))
  e.ab <- ergm(samplike~nodeofactor(function(x) x %v% "Trinity", base=(2:3)), estimate="MPLE")
  expect_summary(s.a, e.a, c(18,36), -c(1.0217, .8353))
  expect_summary(s.ab, e.ab, c(31,30,27), -.8287)
})

test_that("odegree, directed", {
  s.d <- summary(samplike~odegree(2:3))
  e.d <- ergm(samplike~odegree(3), estimate="MPLE")
  s.db <- summary(samplike~odegree(1:3, function(x) x %v% "group"))
  e.db <- ergm(samplike~odegree(4, "group"), estimate="MPLE")
  s.dbh <- summary(samplike~odegree(4:5, ~group, TRUE))
  e.dbh <- ergm(samplike~odegree(2, "group", TRUE), estimate="MPLE")
  expect_summary(s.d, e.d, c(0,1), -.1625189)
  expect_summary(s.db, e.db, c(0,0,0,0,0,1,0,0,0), -c(-1.6292, 0.1112, 0.1625))
  expect_summary(s.dbh, e.dbh, c(6,3), -1.344)
})

test_that("ostar, directed", {
  s.k <- summary(samplike~ostar(1:3))
  e.k <- ergm(samplike~ostar(c(2,4)), estimate="MPLE")
  s.ka <- summary(samplike~ostar(2, "group"))
  e.ka <- ergm(samplike~ostar(2, function(x) x %v% "group"), estimate="MPLE")
  expect_summary(s.k, e.k, c(88,178, 191), c(.1224, -.1986))
  expect_summary(s.ka, e.ka, 88, .1466)
})

test_that("odegree1.5, directed", {
  s.0 <- summary(samplike~odegree1.5)
  e.0 <- ergm(samplike~odegree1.5, estimate="MPLE")
  expect_summary(s.0, e.0, 196.9432, -0.2909)
})

test_that("odegrange, directed", {
  s.0 <- summary(samplike~odegrange(5:8))
  e.0 <- ergm(samplike~odegrange(5:8), estimate="MPLE")
  s.h <- summary(samplike~odegrange(5:8, by=function(x) x %v% "group", homophily=TRUE))
  e.h <- ergm(samplike~odegrange(5:8, by="group", homophily=TRUE), estimate="MPLE")
  expect_summary(s.0, e.0, c(12, 5, 0, 0), -c(0.619, 1.030, Inf, Inf))
  expect_summary(s.h, e.h, c(3, 0, 0, 0), -c(-0.2231, Inf , Inf, Inf ))
})

test_that("receiver, directed", {
  s.0 <- summary(samplike~receiver)
  e.0 <- ergm(samplike~receiver, estimate="MPLE")
  s.b <- summary(samplike~receiver(nodes=-(2:16)))
  e.b <- ergm(samplike~receiver(nodes=-(3:18)), estimate="MPLE")
  expect_summary(s.0, e.0,
    c(8, 4, 2, 5, 3, 5, 7, 11, 10, 6, 3, 6, 3, 5, 3, 2, 3),
    c(-0.1178,-1.1787,-2.0149,-0.8755,-1.5404,-0.8755,
      -0.3567, 0.6061, 0.3567,-0.6061,-1.5404,-0.6061,
      -1.5404,-0.8755,-1.5404,-2.0149,-1.5404))
  expect_summary(s.b, e.b, c(2,2,3), c(-2.0149, -0.1178))
})

test_that("sender, directed", {
  s.0 <- summary(samplike~sender)
  e.0 <- ergm(samplike~sender, estimate="MPLE")
  s.b <- summary(samplike~sender(nodes=-(2:16)))
  e.b <- ergm(samplike~sender(nodes=-(3:18)), estimate="MPLE")
  expect_summary(s.0, e.0,
    c(5, 4, 4, 4, 5, 6, 4, 6, 5, 5, 6, 5, 5, 3, 5, 4, 6),
    -c(0.8755,1.1787,1.1787,1.1787,0.8755,0.6061,1.1787,
       0.6061,0.8755,0.8755,0.6061,0.8755,0.8755,1.5404,
       0.8755,1.1787,0.6061))
  expect_summary(s.b, e.b, c(6,4,6), -c(.6061, .8755))
})

test_that("simmelian, directed", {
  s.0 <- summary(samplike~simmelian)
  e.0 <- ergm(samplike~simmelian, estimate="MPLE")
  expect_summary(s.0, e.0, 8, .6069)
})

test_that("simmelianties, directed", {
  s.0 <- summary(samplike~simmelianties)
  e.0 <- ergm(samplike~simmelianties, estimate="MPLE")
  expect_summary(s.0, e.0, 32, .1984)
})

test_that("transitive, directed", {
  s.0 <- summary(samplike~transitive)
  e.0 <- ergm(samplike~transitive, estimate="MPLE")
  expect_summary(s.0, e.0, 154, -.07745)
})

test_that("transitiveties, directed", {
  s.0 <- summary(samplike~transitiveties)
  e.0 <- ergm(samplike~transitiveties, estimate="MPLE")
  expect_summary(s.0, e.0, 69, -.4116)
})

test_that("ttriple=ttriad, directed", {
  s.0 <- summary(samplike~ttriple)
  e.0 <- ergm(samplike~ttriad, estimate="MPLE")
  s.a <- summary(samplike~ttriple("group"))
  e.a <- ergm(samplike~ttriple(~group), estimate="MPLE")
  s.ad <- summary(samplike~ttriad("group", diff=TRUE))
  e.ad <- ergm(samplike~ttriple(function(x) x %v% "group", diff=TRUE), estimate="MPLE")
  expect_summary(s.0, e.0, 154, -.07745)
  expect_summary(s.a, e.a, 121, .09518)
  expect_summary(s.ad, e.ad, c(26,14,81), c(-.05078, .38935, .13469))
})
