#  File tests/testthat/test-ergm-san.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

context("test-ergm-san.R")

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
