#  File tests/testthat/test-predict.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

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




test_that("predict.formula(net ~ edges + offset(edges))", {
  net <- network.initialize(4, directed=FALSE)

  # edges + offset(edges)
  expect_silent(
    # Odds = 1/4 * 4 = 1
    # P = 0.5
    p <- predict(net ~ edges + offset(edges), c(log(1/4),  log(4)))
  )
  expect_equal(p$p, rep(0.5, 6))
  
  net[1,2:4] <- 1
  expect_equal(
    predict(ergm(net ~ edges + offset(edges), offset.coef=log(4)))$p,
    rep(0.5, 6)
  )
})


test_that("predict.formula(net ~ edges + offset(nodematch))", {
  net <- network.initialize(4, directed=FALSE)
  net %v% "a" <- a <- c(1,1,2,2)

  expect_silent(
    p <- predict(
      net ~ edges + offset(nodematch("a", diff=FALSE)), 
      c(log(1/4),  log(4))
    )
  )
  match_on_a <-  a[p$tail] == a[p$head]
  expect_equal(
    p$p,
    ifelse(match_on_a, 0.5, 0.2)
  )
})




test_that("predict.formula(net ~ edges + degree(1)", {
  net <- network.initialize(3, directed=FALSE)
  net[1,2] <- 1
  expect_silent(
    p <- predict(
      # logodds(1--2) = log(1/4) + log(4)*2
      # odds(1--2) = 16/4 = 4
      # P(1--2) = 4/5
      # logodds(1--3 | 2--3) = log(1/4) + log(4) * 0
      # odds(1--3 | 2--3) = 1/4
      # P(1--3 | 2--3) = 1/5
      net ~ edges + degree(1), 
      c(log(1/4),  log(4))
    )
  )
  expect_equal(
    p$p,
    with(p, ifelse(tail == 1 & head == 2, 4/5, 1/5))
  )
})




test_that("it works for offsets and non-finite offset coefs (and MPLE existence check works)", {
  data("faux.mesa.high")
  expect_warning(fit <- ergm(
    faux.mesa.high ~ edges
    + nodefactor("Grade")
    + nodematch("Grade", diff=T)
    + offset(nodematch("Sex", diff = TRUE, levels = c(1, 2))),
    offset.coef = rep(-Inf, 2)
  ), "^The MPLE does not exist!$")
  expect_silent(
    p <- predict(fit)
  )
  expect_true(
    all(is.finite(p$p))
  )
})


test_that("matrix output of predict() is properly named", {
  data(g4)
  set.seed(666)
  fit <- ergm(g4 ~ edges)
  p.cond <- predict(fit, conditional = TRUE, output = "matrix")
  expect_identical(rownames(p.cond), g4 %v% "vertex.names")
  expect_identical(colnames(p.cond), g4 %v% "vertex.names")
  p.uncond <- predict(fit, conditional = FALSE, output = "matrix", nsim = 2)
  expect_identical(rownames(p.uncond), g4 %v% "vertex.names")
  expect_identical(colnames(p.uncond), g4 %v% "vertex.names")
})
