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

test_that("SAN moves from a sparser network to a denser one with desired triadic attributes", {
	x <- network(100, density = 0.05, directed = FALSE)
	y <- san(x ~ edges + triangles, target.stats = c(600, 300))
	z <- summary(y ~ edges + triangles)

	expect_true(z["edges"] > 580 && z["edges"] < 620)
	expect_true(z["triangle"] > 290 && z["triangle"] < 310)
})


test_that("SAN correctly adjusts inward and outward sums while maintaining edge count", {
	x <- network(100, numedges = 100)
	x %v% 'prop' <- runif(100, 0, 2)
	y <- san(x ~ edges + nodeicov('prop') + nodeocov('prop'), target.stats = c(100, 75, 125))
	z <- summary(y ~ edges + nodeicov('prop') + nodeocov('prop'))

	expect_true(z["edges"] > 95 && z["edges"] < 105)
	expect_true(z["nodeicov.prop"] > 70 && z["nodeicov.prop"] < 80)
	expect_true(z["nodeocov.prop"] > 120 && z["nodeocov.prop"] < 130)
})
