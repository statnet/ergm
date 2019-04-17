
context("test-ergm-san.R")

test_that("SAN moves from a sparser network to a denser one with desired triadic attributes", {
	x <- network(100, density = 0.05, directed = F)
	y <- structure(replicate(5, san(x ~ edges + triangles, target.stats = c(600, 300), only.last = F)), class = "network.list")
	z <- summary(y ~ edges + triangles)

	mean_edges <- mean(z[,"edges"])
	mean_triangles <- mean(z[,"triangle"])

	expect_true(mean_edges > 580 && mean_edges < 620)
	expect_true(mean_triangles > 290 && mean_triangles < 310)
})


test_that("SAN correctly adjusts inward and outward sums while maintaining edge count", {
	x <- network(100, numedges = 100)
	x %v% 'prop' <- runif(100, 0, 2)
	y <- structure(replicate(5, san(x ~ edges + nodeicov('prop') + nodeocov('prop'), target.stats = c(100, 75, 125), only.last = F)), class = "network.list")
	z <- summary(y ~ edges + nodeicov('prop') + nodeocov('prop'))

	mean_edges <- mean(z[,"edges"])
	mean_nodeicov <- mean(z[,"nodeicov.prop"])
	mean_nodeocov <- mean(z[,"nodeocov.prop"])


	expect_true(mean_edges > 95 && mean_edges < 105)
	expect_true(mean_nodeicov > 70 && mean_nodeicov < 80)
	expect_true(mean_nodeocov > 120 && mean_nodeocov < 130)
})