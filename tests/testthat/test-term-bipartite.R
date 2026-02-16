#  File tests/testthat/test-term-bipartite.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

# a bipartite nw
set.seed(143)
b1 <- floor(runif(60, 1,100))
b2 <- floor(runif(60, 101, 130))
exbip.el <- cbind(b1,b2)
bipnw <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw %v% "Letter" <- letters[1:3]
bipnw %v% "Cost" <- c(3,2,1)


# another bipartite nw with more ties and 2 attributes
set.seed(258)
b1 <- floor(runif(150, 1,200))
b2 <- floor(runif(150, 201, 400))
exbip.el <- cbind(b1,b2)
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=200, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color

test_that("b1concurrent, bipartite, undirected", {
  s.0 <- summary(bipnw~b1concurrent)
  e.0 <- ergm(bipnw~b1concurrent, estimate="MPLE")
  s.b <- summary(bipnw~b1concurrent("Letter"))
  e.b <- ergm(bipnw~b1concurrent(function(x) x %v% "Letter"), estimate="MPLE")
  expect_summary(s.0, e.0, 12, -3.961)
  expect_summary(s.b, e.b, c(4,5,3), c(-4.143, -3.62, -4.105))
})

test_that("b1cov, bipartite, undirected", {
  s.a <- summary(bipnw~b1cov("Cost"))
  e.a <- ergm(bipnw~b1cov(~Cost), estimate="MPLE")
  s.at <- summary(bipnw~b1cov(~poly(Cost,2,raw=TRUE)))
  expect_summary(s.a, e.a, 121, -2.212)
  expect_equal(s.at, c(121, 283), ignore_attr=TRUE)
})

test_that("b1degrange, bipartite, undirected", {
  s.d <- summary(bipnw~b1degrange(from=c(0,1,2),to=c(3,3,Inf)))
  e.d <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf)), estimate="MPLE")
  s.anh <- summary(bipnw~b1degrange(from=c(0,1,2),to=c(3,3,Inf),by="Letter",homophily=FALSE))
  e.anh <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=FALSE), estimate="MPLE")
  s.dh <- summary(bipnw~b1degrange(from=c(0,1,2),to=c(3,3,Inf),by="Letter",homophily=TRUE))
  e.dh <- ergm(bipnw~b1degrange(from=c(1,2),to=c(Inf,Inf),by=function(x) x %v% "Letter",homophily=TRUE), estimate="MPLE")
  expect_summary(s.d, e.d, c(96,38,12), c(-4.027, -3.961))
  expect_summary(s.anh, e.anh, c(32,11,4,31,11,5,33,16,3), c(-4.215, -4.143, -4.284, -3.620, -3.636, -4.105))
  expect_summary(s.dh, e.dh, c(100,19,3), c(-3.891, -3.143))
})

test_that("b1degree, bipartite, undirected", {
  s.d <- summary(bipnw~b1degree(0:3))
  e.d <- ergm(bipnw~b1degree(1:3), estimate="MPLE")
  s.db <- summary(bipnw~b1degree(c(0,2:4), by="Letter"))
  e.db <- ergm(bipnw~b1degree(2, by=~Letter), estimate="MPLE")

  expect_summary(s.d, e.d, c(58,30,8,2), -c(2.991, 5.442, 6.484))
  expect_summary(s.db, e.db, c(21,2,1,1,20,3,1,1,17,3,0,0), -c(1.481, .959, 1.431))
})

test_that("b1dsp, bipartite", {
  s.d0 <- summary(bipnw~b1dsp(0))
  e.d0 <- ergm(bipnw~b1dsp(0), estimate="MPLE")
  s.d1 <- summary(bipnw~b1dsp(1:3))
  e.d1 <- ergm(bipnw~b1dsp(1:3), estimate="MPLE")

  expect_summary(s.d0, e.d0, 4900, 2.11343)
  expect_summary(s.d1, e.d1, c(49,1,0), -c(2.096629, 3.0399271))
  expect_true(is.infinite(coef(e.d1)[3]))
})

test_that("b1factor, bipartite, undirected", {
  s.a <- summary(bipnw~b1factor("Letter"))
  e.a <- ergm(bipnw~b1factor(function(x) x %v% "Letter"), estimate="MPLE")
  s.ab <- summary(bipnw~b1factor(~Letter, levels=-3))
  e.ab <- ergm(bipnw~b1factor("Letter", base=2), estimate="MPLE")
  expect_summary(s.a, e.a, c(21,19), -c(3.797, 3.899))
  expect_summary(s.ab, e.ab, c(20,21), -c(3.877, 3.899))
})

test_that("b1mindegree, bipartite, undirected", {
  s.d <- summary(bipnw~b1mindegree(1:3))
  e.d <- ergm(bipnw~b1mindegree(1:3), estimate="MPLE")
  expect_summary(s.d, e.d, c(42,12,4), -c(4.027, 3.961, 3.584))
})

test_that("b1sociality, bipartite, undirected", {
  s.d <- summary(bipnw~b1sociality(nodes=94:96))
  e.d <- ergm(bipnw~b1sociality(nodes=94:96), estimate="MPLE")
  expect_summary(s.d, e.d, c(4,4,2), -c(1.833, 1.833, 2.603))
})

test_that("b1star, bipartite, undirected", {
  s.k <- summary(bipnw~b1star(1:2))
  e.k <- ergm(bipnw~b1star(1:2), estimate="MPLE")
  s.ka <- summary(bipnw~b1star(2:3, function(x) x %v% "Letter"))
  e.ka <- ergm(bipnw~b1star(2:2, ~Letter), estimate="MPLE")
  expect_summary(s.k, e.k, c(60,26), -c(4.0823, -.3179))
  expect_summary(s.ka, e.ka, c(3,0), -3.157)
})

test_that("b1starmix, bipartite, undirected", {
  s.ka <- summary(bipnw2~b1starmix(2, "Letter"))
  e.ka <- ergm(bipnw2~b1starmix(2, "Letter"), estimate="MPLE")
  s.kab <- summary(bipnw2~b1starmix(1, "Letter", base=2))
  e.kab <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3), estimate="MPLE")
  s.kad <- summary(bipnw2~b1starmix(1, "Letter", diff=FALSE))
  e.kad <- ergm(bipnw2~b1starmix(1, "Letter", diff=FALSE), estimate="MPLE")
  s.kabd <- summary(bipnw2~b1starmix(1, "Letter", base=2, diff=FALSE))
  e.kabd <- ergm(bipnw2~b1starmix(1, "Letter", base=2:3, diff=FALSE), estimate="MPLE")
  expect_summary(s.ka, e.ka, c(9,4,7,4), -c(4.870, 5.802, 5.267, 5.962))
  expect_summary(s.kab, e.kab, c(36, 39, 40), -c(5.613, 5.497))
  expect_summary(s.kad, e.kad, c(75, 75), -c(5.567, 5.567))
  expect_summary(s.kabd, e.kabd, 75, -5.567)
})

test_that("b1twostar, bipartite, undirected", {
  s.a <- summary(bipnw2~b1twostar("Letter"))
  e.a <- ergm(bipnw2~b1twostar(function(x) x %v% "Letter"), estimate="MPLE")
  s.aa <- summary(bipnw2~b1twostar("Letter", "Color"))
  e.aa <- ergm(bipnw2~b1twostar(~Letter, "Color"), estimate="MPLE")
  s.ab <- summary(bipnw2~b1twostar(function(x) x %v% "Letter", levels2=-(2:4)))
  e.ab <- ergm(bipnw2~b1twostar("Letter", levels2=-c(1,3,5)), estimate="MPLE")
  s.aab <- summary(bipnw2~b1twostar(~Letter, "Color", base=(2:4)))
  e.aab <- ergm(bipnw2~b1twostar("Letter", "Color", base=c(1,3,5)), estimate="MPLE")
  expect_summary(s.a, e.a, c(9,4,15,17,7,4), -c(4.523, 5.22, 4.773, 4.593, 4.881, 5.446))
  expect_summary(s.aa, e.aa, c(13,2,13,15,5,8), -c(4.548, 6.281, 4.882, 4.702, 4.758, 4.364))
  expect_summary(s.ab, e.ab, c(9,7,4), -c(5.22, 4.593, 5.446))
  expect_summary(s.aab, e.aab, c(13,5,8), -c(6.281, 4.702, 4.364))
})

test_that("b2concurrent, bipartite, undirected", {
  s.0 <- summary(bipnw~b2concurrent)
  e.0 <- ergm(bipnw~b2concurrent, estimate="MPLE")
  s.b <- summary(bipnw~b2concurrent(function(x) x %v% "Letter"))
  e.b <- ergm(bipnw~b2concurrent("Letter"), estimate="MPLE")
  expect_summary(s.0, e.0, 20, -3.497)
  expect_summary(s.b, e.b, c(8,6,6), -c(2.803, 4.190, 2.803))
})

test_that("b2cov, bipartite, undirected", {
  s.a <- summary(bipnw~b2cov("Cost"))
  e.a <- ergm(bipnw~b2cov(~Cost), estimate="MPLE")
  s.at <- summary(bipnw~b2cov(~poly(Cost,2,raw=TRUE)))
  expect_summary(s.a, e.a, c(129), -c(2.191))
  expect_equal(s.at, c(129,317), ignore_attr=TRUE)
})

test_that("b2degrange, bipartite, undirected", {
  s.d <- summary(bipnw~b2degrange(from=c(0,1,2),to=c(3,3,Inf)))
  e.d <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf)), estimate="MPLE")
  s.anh <- summary(bipnw~b2degrange(from=c(0,1,2),to=c(3,3,Inf),by=function(x) x %v% "Letter",homophily=FALSE))
  (e.anh <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=~Letter,homophily=FALSE), estimate="MPLE")) |>
    expect_warning("The MPLE does not exist!")
  s.dh <- summary(bipnw~b2degrange(from=c(0,1,2),to=c(3,3,Inf),by=function(x) x %v% "Letter",homophily=TRUE))
  e.dh <- ergm(bipnw~b2degrange(from=c(1,2),to=c(Inf,Inf),by=~Letter,homophily=TRUE), estimate="MPLE")
  expect_summary(s.d, e.d, c(18,15,20), -c(3.912,3.497))
  expect_summary(s.anh, e.anh, c(4,4,8,7,7,6,7,4,6), -c(-13.566, 2.803, -14.365, 4.190, 5.704, 2.803))
  expect_summary(s.dh, e.dh, c(29,19,3), -c(3.03, 4.46))
})

test_that("b2degree, bipartite, undirected", {
  s.d <- summary(bipnw~b2degree(0:3))
  e.d <- ergm(bipnw~b2degree(1:3), estimate="MPLE")
  s.db <- summary(bipnw~b2degree(0:3, by="Letter"))
  e.db <- ergm(bipnw~b2degree(2, by=~Letter), estimate="MPLE")
  expect_summary(s.d, e.d, c(3,6,9,8), c(1.7203, 1.4941, .6768))
  expect_summary(s.db, e.db, c(0,1,3,2,0,4,3,3,3,1,3,3), c(1.0498, -.3001, 1.0217))
})

test_that("b2dsp, bipartite", {
  s.d0 <- summary(bipnw~b2dsp(0))
  e.d0 <- ergm(bipnw~b2dsp(0), estimate="MPLE")
  s.d1 <- summary(bipnw~b2dsp(1:3))
  e.d1 <- ergm(bipnw~b2dsp(1:3), estimate="MPLE")

  expect_summary(s.d0, e.d0, 381, 2.829767)
  expect_summary(s.d1, e.d1[1:2], c(24,1,0), -c(2.804156, 5.140782))
  expect_true(is.infinite(coef(e.d1)[3]))
})

test_that("b2mindegree, bipartite, undirected", {
  s.d <- summary(bipnw~b2mindegree(1:3))
  e.d <- ergm(bipnw~b2mindegree(1:3), estimate="MPLE")
  expect_summary(s.d, e.d, c(26,20,11), -c(3.912,3.497,3.604))
})

test_that("b2factor, bipartite, undirected", {
  s.a <- summary(bipnw~b2factor("Letter"))
  e.a <- ergm(bipnw~b2factor(function(x) x %v% "Letter"), estimate="MPLE")
  s.ab <- summary(bipnw~b2factor(~Letter, levels=-3))
  e.ab <- ergm(bipnw~b2factor("Letter", base=2), estimate="MPLE")
  expect_summary(s.a, e.a, c(19,16), -c(3.944, 4.119))
  expect_summary(s.ab, e.ab, c(25,19), -c(3.555, 4.119))
})

test_that("b2sociality, bipartite, undirected", {
  s.d <- summary(bipnw~b2sociality(nodes=2:6))
  e.d <- ergm(bipnw~b2sociality(nodes=2:6), estimate="MPLE")
  expect_summary(s.d, e.d, c(2, 3, 2, 3, 1), -c(3.892, 3.476, 3.892, 3.476, 4.595))
})

test_that("b2star, bipartite, undirected", {
  s.k <- summary(bipnw~b2star(1:2))
  e.k <- ergm(bipnw~b2star(1:2), estimate="MPLE")
  s.ka <- summary(bipnw~b2star(2:3, function(x) x %v% "Letter"))
  e.ka <- ergm(bipnw~b2star(2:2, ~Letter), estimate="MPLE")
  expect_summary(s.k, e.k, c(60,51), -c(3.3457, .2724))
  expect_summary(s.ka, e.ka, c(3,0), -4.464)
})

test_that("b2starmix, bipartite, undirected", {
  s.ka <- summary(bipnw2~b2starmix(2, "Letter"))
  e.ka <- ergm(bipnw2~b2starmix(1, "Letter"), estimate="MPLE")
  s.kab <- summary(bipnw2~b2starmix(1, "Letter", base=2))
  e.kab <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3), estimate="MPLE")
  s.kad <- summary(bipnw2~b2starmix(1, "Letter", diff=FALSE))
  e.kad <- ergm(bipnw2~b2starmix(1, "Letter", diff=FALSE), estimate="MPLE")
  s.kabd <- summary(bipnw2~b2starmix(1, "Letter", base=2, diff=FALSE))
  e.kabd <- ergm(bipnw2~b2starmix(1, "Letter", base=2:3, diff=FALSE), estimate="MPLE")
  expect_summary(s.ka, e.ka, c(6, 8, 3, 6), -c(5.613, 5.641, 5.523, 5.497))
  expect_summary(s.kab, e.kab, c(36, 39, 40), -c(5.613, 5.497))
  expect_summary(s.kad, e.kad, c(71,79), -c(5.627, 5.510))
  expect_summary(s.kabd, e.kabd, 71, -5.627)
})

test_that("b2twostar, bipartite, undirected", {
  s.a <- summary(bipnw2~b2twostar("Letter"))
  e.a <- ergm(bipnw2~b2twostar(~Letter), estimate="MPLE")
  s.aa <- summary(bipnw2~b2twostar(function(x) x %v% "Letter", "Color"))
  e.aa <- ergm(bipnw2~b2twostar("Letter", "Color"), estimate="MPLE")
  s.ab <- summary(bipnw2~b2twostar(~Letter, levels2=-(2:4)))
  e.ab <- ergm(bipnw2~b2twostar(function(x) x %v% "Letter", base=c(1,3,5)), estimate="MPLE")
  s.aab <- summary(bipnw2~b2twostar(~Letter, "Color", base=(2:4)))
  e.aab <- ergm(bipnw2~b2twostar("Letter", "Color", levels2=-c(1,3,5)), estimate="MPLE")
  expect_summary(s.a, e.a, c(6,3,16,16,8,6), -c(5, 5.754, 4.603, 4.780, 4.462, 5.055))
  expect_summary(s.aa, e.aa, c(8, 1, 17, 15, 4, 10), -c(4.823, 6.739, 4.690, 4.702, 5.351, 4.37))
  expect_summary(s.ab, e.ab, c(6,8,6), -c(5.754, 4.78, 5.055))
  expect_summary(s.aab, e.aab, c(8,4,10), -c(6.739, 4.702, 4.370))
})

test_that("coincidence, bipartite, undirected, no test on fitting due to the number of coef is large", {
  s.c <- table(summary(bipnw~coincidence,active=0))
  #e.c <- table(round(ergm(bipnw~coincidence, estimate="MPLE")$coef,0))
  expect_true(all(s.c == c(381,24,1)))
})

test_that("gwb1degree, bipartite", {
  s.d <- summary(bipnw~gwb1degree())
  expect_error(summary(bipnw~gwb1degree(cutoff=1)), ".*Term .gwb1degree. has encountered a network for which mode-1 degree of some node exceeded the cut-off setting of 1. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  e.d <- ergm(bipnw~gwb1degree(.4, fixed=TRUE), estimate="MPLE")
  s.df <- summary(bipnw~gwb1degree(.3, fixed=TRUE))
  e.df <- ergm(bipnw~gwb1degree(.2, fixed=TRUE), estimate="MPLE")
  s.da <- summary(bipnw~gwb1degree(.1, fixed=TRUE, attr="Letter"))
  e.da <- ergm(bipnw~gwb1degree(.1, attr=function(x) x %v% "Letter", fixed=TRUE), estimate="MPLE")
  s.dfa <- summary(bipnw~gwb1degree(.1, TRUE, ~Letter))
  e.dfa <- ergm(bipnw~gwb1degree(.1, TRUE, "Letter"), estimate="MPLE")
  expect_equal(s.d, setNames(summary(bipnw~b1degree(1:29)), paste0("gwb1degree#",1:29)))
  expect_equal(coef(e.d), -6.979, tolerance=0.001, ignore_attr=TRUE)
  expect_summary(s.df, e.df, 45.4137, -8.057)
  expect_equal(coef(e.da), -c(6.729, 6.762, 6.418), ignore_attr=TRUE, tolerance=0.001)
  expect_summary(s.dfa, e.dfa, c(13.39962, 13.49479, 16.28549), -c(6.729, 6.762, 6.418))
})

test_that("gwb1dsp, bipartite", {
  s.d0 <- summary(bipnw~gwb1dsp)
  expect_silent(summary(bipnw~gwb1dsp(cutoff=2)))
  expect_error(summary(bipnw~gwb1dsp(cutoff=1)), ".*Term .dgwb1dsp. has encountered a network for which number of outgoing shared partners on some dyad exceeded the cut-off setting of 1. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  expect_warning(summary(bipnw~gwb1dsp(.3)), ".*When .fixed=FALSE. parameter .decay. has no effect..*")
  s.d2 <- summary(bipnw~gwb1dsp(.3, TRUE))
  e.d2 <- ergm(bipnw~gwb1dsp(.3, TRUE), estimate="MPLE")

  expect_equal(s.d0, setNames(c(49,1,rep(0,27)), paste0("b1dsp#", seq_along(s.d0))))
  expect_summary(s.d2, e.d2, c(gwb1dsp.fixed.0.3=50.25918), c(gwb1dsp.fixed.0.3=-2.105815))
})

test_that("gwb2degree, bipartite", {
  s.d <- summary(bipnw~gwb2degree())
  expect_error(summary(bipnw~gwb2degree(cutoff=2)), ".*Term .gwb2degree. has encountered a network for which mode-2 degree of some node exceeded the cut-off setting of 2. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  e.d <- ergm(bipnw~gwb2degree(.4, fixed=TRUE), estimate="MPLE")
  s.df <- summary(bipnw~gwb2degree(.3, fixed=TRUE))
  e.df <- ergm(bipnw~gwb2degree(.2, fixed=TRUE), estimate="MPLE")
  s.da <- summary(bipnw~gwb2degree(.1, fixed=TRUE, attr="Letter"))
  e.da <- ergm(bipnw~gwb2degree(.1, attr=function(x) x %v% "Letter", fixed=TRUE), estimate="MPLE")
  s.dfa <- summary(bipnw~gwb2degree(.1, TRUE, ~Letter))
  e.dfa <- ergm(bipnw~gwb2degree(.1, TRUE, "Letter"), estimate="MPLE")
  expect_equal(s.d, setNames(summary(bipnw~b2degree(1:30)), paste0("gwb2degree#",1:30)))
  expect_equal(coef(e.d), -25.99385, tolerance=0.001, ignore_attr=TRUE)
  expect_summary(s.df, e.df, 31.97479, -32.78813)
  expect_equal(coef(e.da), -c(33.82191, 24.76756, 34.28992), ignore_attr=TRUE, tolerance=0.001)
  expect_summary(s.dfa, e.dfa, c(9.809166, 10.598143, 7.598143), -c(33.82191, 24.76756, 34.28992))
})

test_that("gwb2dsp, bipartite", {
  s.d0 <- summary(bipnw~gwb2dsp)
  expect_silent(summary(bipnw~gwb2dsp(cutoff=2)))
  expect_error(summary(bipnw~gwb2dsp(cutoff=1)), ".*Term .dgwb2dsp. has encountered a network for which number of incoming shared partners on some dyad exceeded the cut-off setting of 1. This can usually be remedied by increasing the value of the term argument .cutoff. or the corresponding term option .gw.cutoff...*")
  expect_warning(summary(bipnw~gwb2dsp(.3)), ".*When .fixed=FALSE. parameter .decay. has no effect..*")
  s.d2 <- summary(bipnw~gwb2dsp(.3, TRUE))
  e.d2 <- ergm(bipnw~gwb2dsp(.3, TRUE), estimate="MPLE")

  expect_equal(s.d0, setNames(c(24,1,rep(0,28)), paste0("b2dsp#", seq_along(s.d0))))
  expect_summary(s.d2, e.d2, c(gwb2dsp.fixed.0.3=25.25918), c(gwb2dsp.fixed.0.3=-2.875923))
})
