#  File tests/testthat/test-ergm-san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################


n <- 50

test_that("SAN moves from a sparser network to a denser one with desired triadic attributes", {
	x <- network(n, density = 0.05/100*n, directed = FALSE)
	y <- san(x ~ edges + triangles, target.stats = c(n*6, n*3))
	z <- summary(y ~ edges + triangles)

	expect_true(z["edges"] > n*5.8 && z["edges"] < n*6.2)
	expect_true(z["triangle"] > n*2.9 && z["triangle"] < n*3.1)
})


test_that("SAN correctly adjusts inward and outward sums while maintaining edge count", {
	x <- network(n, numedges = n)
	x %v% 'prop' <- runif(n, 0, 2)
	y <- san(x ~ edges + nodeicov('prop') + nodeocov('prop'), target.stats = c(n, n*.75, n*1.25))
	z <- summary(y ~ edges + nodeicov('prop') + nodeocov('prop'))

	expect_true(z["edges"] > n*.95 && z["edges"] < n*1.05)
	expect_true(z["nodeicov.prop"] > n*.7 && z["nodeicov.prop"] < n*.8)
	expect_true(z["nodeocov.prop"] > n*1.2 && z["nodeocov.prop"] < n*1.3)
})

test_that("SAN matches target stats while respecting infinite offsets", {
    x <- network(n, directed=FALSE,density=0)
    x %v% "sex" <- sample(c("M","F"),n,rep=TRUE)
    y <- san(x ~ edges + offset(nodematch("sex")), target.stats=c(6*n), offset.coef=c(-Inf))
    z <- summary(y ~ edges + offset(nodematch("sex")))
    expect_true(z["edges"] > 5.9*n && z["edges"] < 6.1*n)
    expect_true(z["offset(nodematch.sex)"] == 0)
})

test_that("SAN matches target stats while respecting infinite dyad-dependent offsets", {
    x <- network(n, directed=FALSE,density=0)
    y <- san(x ~ edges + offset(concurrent), target.stats=c(floor(n/2)), offset.coef=c(-Inf))
    z <- summary(y ~ edges + offset(concurrent))
    expect_true(z["edges"] >= 0.95*floor(n/2))
    expect_true(z["offset(concurrent)"] == 0)
})

test_that("weighted SAN matches target stats while respecting infinite offsets", {
    x <- network(n, directed=FALSE, numedges=0)
    x %v% "sex" <- rep(c("M","F"),length.out=n)
    y <- san(x ~ sum + offset(nodematch("sex")), reference=~Unif(0,n/5),target.stats=c(n^3/160),response="ea",offset.coef=c(-Inf))
    z <- summary(y ~ sum + offset(nodematch("sex")), response="ea")
    expect_true(z["sum"] > .98*n^3/160 && z["sum"] < 1.02*n^3/160)
    expect_true(z["offset(nodematch.sum.sex)"] == 0)
})

test_that("SAN errors when passed the wrong number of offsets", {
    x <- network(n, directed=FALSE,density=0)
    expect_error(san(x ~ edges + offset(concurrent), target.stats=c(6*n)), paste0("Length of ", sQuote("offset.coef"), " in SAN is 0, while the number of offset coefficients in the model is 1."))
    expect_error(san(x ~ edges + offset(concurrent), target.stats=c(6*n), offset.coef=c(1,2)), paste0("Length of ", sQuote("offset.coef"), " in SAN is 2, while the number of offset coefficients in the model is 1."))
})

test_that("san.ergm does not default to offsets in the ergm", {
    x <- network(n, directed=FALSE,numedges=1)
    e <- ergm(x ~ edges + offset(concurrent), offset.coef=c(-Inf),estimate="MPLE")
    expect_error(san(e, target.stats=c(floor(n/2))))
    y <- san(e, target.stats=c(floor(n/2)), offset.coef=c(-Inf))
    z <- summary(y ~ edges + offset(concurrent))
    expect_true(z["edges"] >= 0.95*floor(n/2))
    expect_true(z["offset(concurrent)"] == 0)
})

test_that("SAN works with curved terms", {
    x <- network(n, directed=FALSE,numedges=1)
    y <- san(x ~ edges + gwesp(0,fixed=T), target.stats=c(100,10))
    z <- summary(y ~ edges + gwesp(0,fixed=T))
    expect_true(z["edges"] >= 98 && z["edges"] <= 102)
    expect_true(z["gwesp.fixed.0"] >= 9 && z["gwesp.fixed.0"] <= 11)
    
    e <- ergm(x ~ edges + offset(gwesp(0,fixed=T)), offset.coef=c(-Inf), estimate="MPLE")
    y <- san(e, target.stats=c(250), offset.coef=c(-Inf))
    z <- summary(y ~ edges + offset(gwesp(0,fixed=T)))
    expect_true(z["edges"] >= 245 && z["edges"] <= 255)
    expect_true(z["offset(gwesp.fixed.0)"] == 0)
    
    y <- san(x ~ edges + gwesp(cutoff=2), target.stats=c(500,20,10))
    z <- summary(y ~ edges + gwesp(cutoff=2))
    expect_true(z["edges"] >= 495 && z["edges"] <= 505)
    expect_true(z["esp#1"] >= 19 && z["esp#1"] <= 21)
    expect_true(z["esp#2"] >= 9 && z["esp#2"] <= 11)
    
    e <- ergm(x ~ edges + offset(degree(3)) + gwesp(0,fixed=T), offset.coef=c(-Inf), estimate="MPLE")
    y <- san(e, target.stats=c(30,9), offset.coef=c(-Inf))
    z <- summary(y ~ edges + gwesp(0,fixed=T))
    expect_true(z["edges"] >= 29 && z["edges"] <= 31)
    expect_true(z["gwesp.fixed.0"] >= 8 && z["gwesp.fixed.0"] <= 10)
    
    e <- ergm(x ~ edges + offset(degree(3)) + gwesp(cutoff=2), offset.coef=c(-Inf), control=control.ergm(MCMLE.maxit=1, loglik.control=control.logLik.ergm(nsteps=1)))
    y <- san(e, target.stats=c(30,9,0), offset.coef=c(-Inf))
    z <- summary(y ~ edges + gwesp(cutoff=2))
    expect_true(z["edges"] >= 29 && z["edges"] <= 31)
    expect_true(z["esp#1"] >= 8 && z["esp#1"] <= 10)
    expect_true(z["esp#2"] == 0)
})

test_that("SAN offsets work with curved terms", {
    x <- network(n, directed=FALSE,numedges=1)
    expect_error(san(x ~ edges + offset(gwesp(1)), target.stats=c(100)))
    expect_error(san(x ~ edges + offset(gwesp(1)), target.stats=c(100), offset.coef=c(1,1)))    
    y <- san(x ~ edges + offset(gwesp(1)), target.stats=c(100), offset.coef=c(1))
    z <- summary(y ~ edges)
    expect_true(z["edges"] >= 98 && z["edges"] <= 102)

    expect_error(san(x ~ edges + offset(gwesp(1)) + triangle, target.stats=c(100, 10)))
    expect_error(san(x ~ edges + offset(gwesp(1)) + triangle, target.stats=c(100, 10), offset.coef=c(1,1)))    
    y <- san(x ~ edges + offset(gwesp(1)) + triangle, target.stats=c(100, 10), offset.coef=c(-Inf))
    z <- summary(y ~ triangle)
    expect_true(z["triangle"] == 0)

    expect_error(san(x ~ edges + gwesp(0.1) + offset(gwnsp(1)), target.stats=c(100)))
    expect_error(san(x ~ edges + gwesp(0.1) + offset(gwnsp(1)), target.stats=c(100), offset.coef=c(1,1)))    
    y <- san(x ~ edges + gwesp(0.1) + offset(gwnsp(1)), target.stats=c(100, 15:1, rep(0, 15)), offset.coef=c(-Inf))
    expect_true(all(summary(y ~ gwnsp(1)) == rep(0, 30)))

    expect_error(san(x ~ edges + offset(gwnsp(1)) + gwesp(0.1), target.stats=c(100)))
    expect_error(san(x ~ edges + offset(gwnsp(1)) + gwesp(0.1), target.stats=c(100), offset.coef=c(1,1)))    
    y <- san(x ~ edges + offset(gwnsp(1)) + gwesp(0.1), target.stats=c(100, 15:1, rep(0, 15)), offset.coef=c(-Inf))
    expect_true(all(summary(y ~ gwnsp(1)) == rep(0, 30)))
})
