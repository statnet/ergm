#  File tests/testthat/test-term-flexible.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

expect_summary <- function(s, e, value, coefficients, tolerance=0.001) {
  expect_equal(s, value, tolerance=tolerance, ignore_attr=TRUE)
  expect_equal(coef(e)[1:length(coefficients)], coefficients, tolerance=tolerance, ignore_attr=TRUE)
}

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
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color


# a directed nw
load("sampson.wrong.RData") # Old (wrong) version of sampson's monks
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")


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

test_that("absdiff, no required type, independent", {
  s.a <- summary(fmh ~ absdiff("Grade"))
  e.a <- ergm(fmh ~ absdiff(function(x) x %v% "Grade"))
  s.ap <- summary(fmh ~ absdiff(~Grade, pow=2))
  e.ap <- ergm(fmh ~ absdiff("Grade", pow=2))

  expect_summary(s.a, e.a, 79, -4.354)
  expect_summary(s.ap, e.ap, 195, -3.41)
})

test_that("absdiffcat, no required type, independent", {
  s.a <- summary(fmh ~ absdiffcat("Grade"))
  e.a <- ergm(fmh ~ absdiffcat("Grade"))
  s.ab <- summary(fmh ~ absdiffcat(function(x) x %v% "Grade", levels=-(4:5)))
  e.ab <- ergm(fmh ~ absdiffcat(~Grade, base=(4:5)))
  expect_summary(s.a, e.a, c(15,15,7,2,1), -c(6.005,5.788,6.063,6.891,6.611))
  expect_summary(s.ab, e.ab, c(15,15,7), -c(6.005,5.788,6.063))
})

test_that("balance, dir or undir", {
  s.0 <- summary(fmh~balance)
  e.0 <- ergm(fmh~balance, estimate="MPLE")
  expect_summary(s.0, e.0, 40139, -.02376)
})

test_that("cycle, either", {
  s.0 <- summary(samplike ~ cycle(2:6))
  e.0 <- ergm(samplike ~ cycle(2:6), estimate="MPLE")
  s.1 <- summary(samplike ~ cycle(3:7,semi=TRUE))
  e.1 <- ergm(samplike ~ cycle(3:7,semi=TRUE), estimate="MPLE")
  s.k <- summary(fmh~cycle(3:6))
  e.k <- ergm(fmh~cycle(c(4,6)), estimate="MPLE")
  expect_summary(s.0, e.0, c(28, 39, 111, 260, 651), c(2.118, -0.539, 0.410, -0.022, -0.049))
  expect_summary(s.1, e.1, c(57, 216, 787, 2908, 10508), c(-0.0091, 0.1439, 0.0704, -0.0311, 0.0011))
  expect_summary(s.k, e.k, c(62,80,138,270), -c(-.1615, .2083))
})

test_that("density, either", {
  s.0 <- summary(fmh~density)
  e.0 <- ergm(samplike~density, estimate="MPLE")
  expect_summary(s.0, e.0, .009708274, -277.5904)
})

test_that("diff, no required type but primarily directed, independent", {
  # Auxiliary variables, useful for calculating the true values of statistics.
  sthd <- outer(samplike%v%"YearsServed",samplike%v%"YearsServed","-") # YS[t]-YS[h]
  sm <- as.matrix(samplike)
  s.a <- summary(samplike ~ diff("YearsServed"))
  e.a <- ergm(samplike ~ diff(~YearsServed))
  s.ad <- summary(samplike ~ diff("YearsServed", dir="h-t"))
  e.ad <- ergm(samplike ~ diff(function(x) x %v% "YearsServed", dir="h-t"))
  s.ads2 <- summary(samplike ~ diff(~YearsServed, sign.action="abs"))
  e.ads2 <- ergm(samplike ~ diff("YearsServed", sign.action="abs"))
  s.ads3 <- summary(samplike ~ diff(~YearsServed, sign.action="pos"))
  e.ads3 <- ergm(samplike ~ diff("YearsServed", sign.action="pos"))
  s.ads4 <- summary(samplike ~ diff(function(x) x %v% "YearsServed", sign.action="neg"))
  e.ads4 <- ergm(samplike ~ diff("YearsServed", sign.action="neg"))
  s.ap <- summary(samplike ~ diff(function(x) x %v% "YearsServed", pow=3))
  e.ap <- ergm(samplike ~ diff("YearsServed", pow=3))

  expect_summary(s.a, e.a, sum(sthd*sm), 0.0631)
  expect_summary(s.ad, e.ad, sum(-sthd*sm), -0.0631)
  expect_summary(s.ads2, e.ads2, sum(abs(sthd)*sm), -0.381)
  expect_summary(s.ads3, e.ads3, sum((sthd+abs(sthd))*sm/2), -0.2843)
  expect_summary(s.ads4, e.ads4, sum((sthd-abs(sthd))*sm/2), 0.504)
  expect_summary(s.ap, e.ap, sum(sthd^3*sm), 0.001844)
})

test_that("dyadcov, either", {
  set.seed(120)
  cov <- matrix(rbinom(324, 1, .5),18,18)
  cov <- cov+t(cov)
  s.x <- summary(samplike~dyadcov(cov))
  e.x <- ergm(samplike ~ dyadcov(cov))
  s.xa <- summary(fmh~dyadcov(fmh, "GradeMet"))
  e.xa <- ergm(fmh ~ dyadcov(fmh, "GradeMet"))
  expect_summary(s.x, e.x, c(31,21,14), -+c(.8546, 1.0732, 1.3467))
  expect_summary(s.xa, e.xa, 641, 12.31787)
})

test_that("edgecov, either", {
  set.seed(64)
  cov <- matrix(rbinom(324, 3, .5),18,18)
  s.x <- summary(samplike~edgecov(cov))
  e.x <- ergm(samplike ~ edgecov(cov))
  s.xa <- summary(samplike~edgecov(samplike, "YearsTrusted"))
  e.xa <- ergm(samplike ~ edgecov(samplike, "YearsTrusted"))
  n.x <- try(summary(samplike~edgecov('dummy')),silent=TRUE)
  set.network.attribute(samplike,'dummy',cov)
  n2.x <- summary(samplike~edgecov('dummy'))
  expect_summary(s.x, e.x, 134, -.5022)
  expect_summary(s.xa, e.xa, 183, Inf)
  expect_true(is(n.x, 'try-error'))
  expect_equal(n2.x, 134, ignore_attr=TRUE)
})

test_that("edges, either", {
  s.0 <- summary(fmh~edges)
  e.0 <- ergm(samplike~edges, estimate="MPLE")
  expect_summary(s.0, e.0, 203, -.9072)
})

run.sp.tests <- function(cache.sp) {
  test_that(paste0("dsp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.d <- summary(fmh~dsp(2:3), cache.sp=cache.sp)
    e.d <- ergm(samplike~dsp(4), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    expect_summary(s.d, e.d, c(75, 23), -.04275)
  })

  test_that(paste0("esp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.d <- summary(fmh~esp(2:3), cache.sp=cache.sp)
    e.d <- ergm(samplike~esp(4), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    expect_summary(s.d, e.d, c(36,13), .3093)
  })

  test_that(paste0("nsp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.d <- summary(fmh~nsp(2:3), cache.sp=cache.sp)
    e.d <- ergm(samplike~nsp(4), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    expect_summary(s.d, e.d, c(39, 10), -+ 1.1096)
  })

  test_that(paste0("gwdsp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.0 <- summary(fmh~gwdsp, cache.sp=cache.sp)
    e.0 <- ergm(samplike~gwdsp(0, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    e.a <- ergm(samplike~gwdsp(.8, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    s.f <- summary(fmh~gwdsp(0, fixed=TRUE), cache.sp=cache.sp)
    s.af <- summary(fmh~gwdsp(.3, fixed=TRUE), cache.sp=cache.sp)
    e.af <- ergm(samplike~gwdsp(.2, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    expect_summary(head(s.0), e.0, c(431, 75, 23, 1, 1, 0), -.3309974)
    expect_equal(coef(e.a), -.1875983, tolerance=0.001, ignore_attr=TRUE)
    expect_equal(s.f, 531, ignore_attr=TRUE)
    expect_summary(s.af, e.af, 558.6369, -.2829672)
  })

  test_that(paste0("gwesp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.0 <- summary(fmh~gwesp, cache.sp=cache.sp)
    e.0 <- ergm(samplike~gwesp(0, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    e.a <- ergm(samplike~gwesp(.8, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    s.f <- summary(fmh~gwesp(0, fixed=TRUE), cache.sp=cache.sp)
    s.af <- summary(fmh~gwesp(.3, fixed=TRUE), cache.sp=cache.sp)
    e.af <- ergm(samplike~gwesp(.2, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))

    expect_summary(head(s.0), e.0, c(70,36,13,0,1,0), -.4115515)
    expect_equal(coef(e.a), -.1898684, tolerance=0.001, ignore_attr=TRUE)
    expect_equal(s.f, 120, ignore_attr=TRUE)
    expect_summary(s.af, e.af, 133.9215, -.3371385)
  })

  test_that(paste0("gwnsp, either, shared partner cache ", if(cache.sp) "enabled" else "disabled"), {
    s.0 <- summary(fmh~gwnsp, cache.sp=cache.sp)
    e.0 <- ergm(samplike~gwnsp(0, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    e.a <- ergm(samplike~gwnsp(.8, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))
    s.f <- summary(fmh~gwnsp(0, fixed=TRUE), cache.sp=cache.sp)
    s.af <- summary(fmh~gwnsp(.3, fixed=TRUE), cache.sp=cache.sp)
    e.af <- ergm(samplike~gwnsp(.2, fixed=TRUE), estimate="MPLE", control=control.ergm(term.options=list(cache.sp=cache.sp)))

    expect_summary(head(s.0), e.0, c(361,39,10,1,0,0), -.4189)
    expect_equal(coef(e.a), -.3123, tolerance=0.001, ignore_attr=TRUE)
    expect_equal(s.f, 411, ignore_attr=TRUE)
    expect_summary(s.af, e.af, 424.7154, -.3934841)
  })
}

run.sp.tests(TRUE)
run.sp.tests(FALSE)

#test_that("hamming, any", {
#  mat.d <- matrix(0,18,18)
#  mat.u <- matrix(0, 205, 205)
#  set.seed(456)
#  # Using a covariate matrix that matches the edges exactly is too easy.
#  cov.d <- cbind(as.edgelist(samplike)[,2:1], rbinom(88, 3, .5))
#  set.seed(145)
#  cov.u <- cbind(as.edgelist(fmh), rbinom(203, 3, .5))
#
#  # although there are 4 non-required inputs, giving
#  # 16 combinations of inputs, I've exlcuded most that
#  # don't involve 'x' because w/o 'x', the results are
#  # 0 or largely negative, as the hamming distance is
#  # compared between identical networks
#  s.0 <- 0# COMMENTED OUT FOR NOW BECAUSE IT'S BROKEN:  summary(samplike~hamming)
#  s.x <- summary(samplike~hamming(mat.d))
#  # and everything commented below is broke.
#
#  # should this really be NA
#  #e.x <- ergm(fmh~hamming(mat.u), estimate="MPLE")
#  ## OK
#  s.xc <- summary(samplike~hamming(mat.d, cov=cov.d))
#  # NA
#  #e.xc <- ergm(fmh~hamming(mat.u, cov=cov.u), estimate="MPLE")
#  # OK
#  s.xd <- summary(samplike~hamming(mat.d, defaultweight=.3))
#  # NA value
#  #e.xd <- ergm(samplike~hamming(mat.d, defaultweight=.3), estimate="MPLE")
#  # OK
#  s.xca <- summary(samplike~hamming(mat.d, cov=samplike, attrname="YearsTrusted"))
#  # NA
#  #e.xca <- ergm(fmh~hamming(mat.u, cov=fmh, attrname="Grade"), estimate="MPLE")
#  # OK
#  s.xcd<- summary(samplike~hamming(mat.d, cov=cov.d, defaultweight=.5))
#  # NA
#  #e.xcd<- ergm(samplike~hamming(mat.d, cov=cov.d, defaultweight=.5), estimate="MPLE")
#  # 0 & NA
#  #s.xcad<- summary(samplike~hamming(mat.d, samplike, "YearsTrusted", .5))
#  #e.xcad<- ergm(samplike~hamming(mat.d, samplike, "YearsTrusted", .5), estimate="MPLE")
#
#  #expect_equal(c(s.0, s.x, s.xc, s.xd, s.xca, s.xcd), c(0, 88, 84, 26.4, 183, 100), ignore_attr=TRUE)
#})

test_that("isolatededges, undirected", {
  s.0 <- summary(fmh~isolatededges)
  e.0 <- ergm(fmh~isolatededges, estimate="MPLE")
  s.1 <- summary(bipnw2~isolatededges)
  e.1 <- ergm(bipnw2~isolatededges, estimate="MPLE")

  expect_summary(s.0, e.0, 4, .01034)
  expect_summary(s.1, e.1, 25, -.1611)
})

test_that("isolates, either", {
  s.0 <- summary(samplike~isolates)
  e.0 <- ergm(fmh~isolates, estimate="MPLE")

  expect_summary(s.0, e.0, 0, 5.10979)
})

#test_that("localtriangle, either", {
#  set.seed(85)
#  x <- matrix(rbinom(324, 2, .5),18,18)
#  s.x <- summary(samplike~localtriangle(x))
#  e.x <- ergm(samplike~localtriangle(x), estimate="MPLE")
#  s.xa <- summary(fmh~localtriangle(fmh, "GradeMet"))
#  expect_summary(s.x, e.x, 56, -.1553)
#  expect_equal(s.xa, 61)
#})

test_that("meandeg, either", {
  s.0 <- summary(samplike~meandeg)
  e.0 <- ergm(fmh~meandeg, estimate="MPLE")
  expect_summary(s.0, e.0, 4.8889, -474.0647)
})

test_that("nodecov, either", {
  s.a <- summary(samplike~nodecov("YearsServed"))
  e.a <- ergm(fmh~nodecov("Grade"), estimate="MPLE")
  s.at <- summary(samplike~nodecov(~YearsServed^2))
  e.at <- ergm(fmh~nodecov(~(.%v%"Grade")^2), estimate="MPLE")
  s.att <- summary(samplike~nodecov(function(x)(x%v%"YearsServed")^2))
  s.attt <- summary(samplike~nodecov(~poly(YearsServed,2,raw=TRUE)))
  expect_summary(s.a, e.a, 906, -.271)
  expect_summary(s.at, e.at, 5036, -.03199)
  expect_equal(s.att, 5036, ignore_attr=TRUE)
  expect_equal(s.attt, c(906, 5036), ignore_attr=TRUE)
})

test_that("nodefactor, either", {
  s.a <- summary(fmh~nodefactor("Grade"))
  e.a <- ergm(samplike~nodefactor(~group), estimate="MPLE")
  s.ab <- summary(fmh~nodefactor(function(x) x %v% "Sex", base=(4:5)))
  e.ab <- ergm(samplike~nodefactor("Trinity", levels=TRUE), estimate="MPLE")
  expect_summary(s.a, e.a, c(75, 65, 36, 49, 28), -c(.9480, .3273))
  expect_summary(s.ab, e.ab, c(235, 171), -c(.4451, .4451, .4706))
})

test_that("nodematch, either", {
  s.a <- summary(fmh~nodematch("Race"))
  e.a <- ergm(samplike~nodematch("Trinity"), estimate="MPLE")
  s.ad <- summary(samplike~nodematch(function(x) x %v% "group", diff=TRUE))
  e.ad <- ergm(fmh~nodematch("Sex", diff=TRUE), estimate="MPLE")
  s.ak <- summary(fmh~nodematch(~Grade, levels=3:4))
  e.ak <- ergm(samplike~nodematch(function(x) x %v% "group", levels=2), estimate="MPLE")
  s.adk <- summary(samplike~nodematch(~Trinity, TRUE, 1:2))
  e.adk <- ergm(fmh~nodematch("Race", TRUE, 2), estimate="MPLE")
  expect_summary(s.a, e.a, 103, -1.45725)
  expect_summary(s.ad, e.ad, c(23,10,30), -c(4.06317, 4.7032))
  expect_summary(s.ak, e.ak, 32, c(1.609, NA))
  expect_summary(s.adk, e.adk, c(8,4), -4.700995)
})

test_that("nodemix, any", {
  s.a <- summary(fmh ~ nodemix("Grade"))
  e.a <- ergm(samplike ~ nodemix(function(x) x %v% "group"), estimate="MPLE")
  s.ab <- summary(bipnw ~ nodemix("Letter", levels2=TRUE))
  e.ab <- ergm(bipnw ~ nodemix(function(x) x %v% "Letter", levels2=-(2:6)))
  s.ab2 <- summary(fmh ~ nodemix("Race", base=1))
  e.ab2 <- ergm(samplike ~ nodemix(~Trinity, base=(3:9)))

  expect_summary(s.a, e.a,
    c(0, 33, 0, 2, 23, 1, 4, 7, 9, 1, 2, 6, 1, 17, 1, 1, 4, 5, 5, 6),
    c(-3.2958369, -2.1747517, -2.5649494, 1.6094379, -3.2958369,  -1.4916549, -1.0986123, 0.9162907))
  expect_summary(s.ab, e.ab, c(9,8,8,7,7,5,4,6,6), -c(3.497, 4.431, 3.989, 3.989))
  expect_summary(s.ab2, e.ab2, c(8,53,13,41,46,0,1,0,0,5,22,10,0,4), -c(1.0116, .82098))
})

test_that("smalldiff", {
  s.ac.d <- summary(samplike~smalldiff("YearsServed", 3))
  s.ac.u <- summary(fmh~smalldiff("Grade", 2))
  s.ac.b <- summary(bipnw~smalldiff("Cost", 1))
  e.ac.d <- ergm(samplike~smalldiff(~YearsServed, 3), estimate="MPLE")
  e.ac.u <- ergm(fmh~smalldiff(~Grade, 2), estimate="MPLE")
  e.ac.b <- ergm(bipnw~smalldiff(function(x) x %v% "Cost", 1), estimate="MPLE")

  expect_summary(s.ac.d, e.ac.d, 78, -.86903)
  expect_summary(s.ac.u, e.ac.u, 193, -4.3525)
  expect_summary(s.ac.b, e.ac.b, 48, -3.8318)
})

test_that("threetrail, either", {
  s.0 <- summary(samplike~threetrail)
  e.0 <- ergm(fmh~threetrail, estimate="MPLE")
  s.k <- summary(samplike~threetrail(levels=2))
  e.k <- ergm(samplike~threetrail(keep=1:2), estimate="MPLE")

  expect_summary(s.0, e.0, c(2103, 2326, 1749, 1897), -.2842)
  expect_summary(s.k, e.k, 2326, -c(.01881, -.00776))
})

test_that("triangles, either", {
  s.0 <- summary(fmh~triangles)
  e.0 <- ergm(samplike~triangles, estimate="MPLE")
  s.a <- summary(fmh~triangles(function(x) x %v% "Race"))
  e.a <- ergm(samplike~triangle("group"), estimate="MPLE")
  s.ad <- summary(samplike~triangles(~Trinity, diff=TRUE))
  e.ad <- ergm(fmh~triangle("Sex", diff=TRUE), estimate="MPLE")

  expect_summary(s.0, e.0, 62, -.06997)
  expect_summary(s.a, e.a, 18, .06354)
  expect_summary(s.ad, e.ad, c(2, 0, 0), -c(.70278, .44099))
})

test_that("triadcensus, either", {
  s.0 <- summary(samplike~triadcensus)
  e.0 <- ergm(fmh~triadcensus, estimate="MPLE")
  s.d <- summary(samplike~triadcensus(3))
  e.d <- ergm(fmh~triadcensus(2:3), estimate="MPLE")

  expect_summary(s.0, e.0, c(205, 190, 12, 24, 24, 68, 34, 5, 0, 35, 15, 6, 5, 18, 8),
    -c(.02559, .06254, -2.61531))
  expect_summary(s.d, e.d, 12, c(-1.749635, 2.228183))
})

test_that("twopath, either", {
  s.0 <- summary(samplike~twopath)
  e.0 <- ergm(fmh~twopath, estimate="MPLE")
  expect_summary(s.0, e.0, 378, -1.297362)
})
