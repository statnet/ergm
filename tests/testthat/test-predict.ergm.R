context("Test predict.formula() and predict.ergm()")

library(ergm)

test_that("predict.formula(type=) give correct results", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(net ~ edges, theta = log(1/5), type="response") # predict.formula()
  )
  expect_silent(
    p.link <- predict(net ~ edges, theta = log(1/5), type="link") # predict.formula()
  )
  expect_true(
    all.equal(p.link$p, log(p.prob$p / (1 - p.prob$p)))
  )
})


test_that("predict.formula(conditional=FALSE)", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(
      net ~ edges, 
      theta = log(1/5),
      nsim = 5,
      type="response", 
      conditional=FALSE
    )
  )
  
})


test_that("works for edges model on small digraph", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    r.f <- predict(net ~ edges, log(1/5)) # predict.formula()
  )
  fit <- ergm(net ~ edges)
  expect_silent(
    r.e <- predict(fit) # predict.ergm()
  )
  expect_true( all.equal(unique(r.f$p), 1/6) )
  expect_identical(
    names(r.f),
    c("tail", "head", "p")
  )
  expect_identical(
    names(r.e),
    c("tail", "head", "p")
  )
  expect_true( all.equal(unique(r.e$p), 1/6) )
})






test_that("predict.formula(output='matrix') works correctly", {
  net <- network.initialize(3, directed=TRUE)
  net[1,2] <- 1
  expect_silent(
    p.prob <- predict(net ~ edges, theta = log(1/5), type="response", output="matrix")
  )
  
})







test_that("works for edges model on small graph", {
  net <- network.initialize(3, directed=FALSE)
  net[1,2] <- 1
  expect_silent(
    r.f <- predict(net ~ edges, log(1/2)) # predict.formula()
  )
  fit <- ergm(net ~ edges)
  expect_silent(
    r.e <- predict(fit) # predict.ergm()
  )
  expect_identical(
    names(r.f),
    c("tail", "head", "p")
  )
  expect_identical(
    names(r.e),
    c("tail", "head", "p")
  )
  expect_true( all.equal(unique(r.f$p), 1/3) )
  expect_true( all.equal(unique(r.f$p), 1/3) )
})




test_that("predict.ergm() works for models with offset terms", {
  net <- network.initialize(4, directed=FALSE)
  net %v% "a" <- c(1,1,2,2)
  expect_silent(
    # Odds = 1/4 * 4 = 1
    # P = 0.5
    p <- predict(net ~ edges + offset(edges), c(log(1/4),  log(4)))
  )
  expect_equal(p$p, rep(0.5, 6))
  # Odds = 1/4 * 2 = 1/2
  # P = 1/3
  expect_silent(
    p <- predict(net ~ edges + offset(edges), c(log(1/4),  log(2)))
  )
  expect_equal(p$p, rep(1/3, 6))
  expect_silent(
    p <- predict(net ~ edges + offset(nodematch("a", diff=FALSE)), c(log(1/4),  log(4)))
  )
  expect_equal(
    p$p,
    c(0.2, 0.5)[c(1,2,1,2,1,1)]
  )
})
