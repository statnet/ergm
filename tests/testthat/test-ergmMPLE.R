#  File tests/testthat/test-ergmMPLE.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
data(faux.mesa.high)
formula <- faux.mesa.high ~ nodematch("Sex")
mplesetup <- ergmMPLE(formula)

test_that("output = 'matrix' works for undirected networks", {
  ord <- order(mplesetup$weights)
  m <- cbind(mplesetup$weights, mplesetup$response, mplesetup$predictor)[ord,]
  expect_equal(m, matrix(c(71, 132, 10284, 10423, 1, 1, 0, 0, 0, 1, 1, 0), 4,3), ignore_attr=TRUE)
})

test_that("output = 'fit' works for undirected networks", {
  modelfit <- ergmMPLE(formula, output="fit")
  alt <- glm(mplesetup$response ~ mplesetup$predictor - 1,
             weights = mplesetup$weights, family="binomial")
  expect_equal(coef(modelfit), coef(alt), ignore_attr=TRUE, tolerance=0.001)
})

data(florentine)

test_that("output = 'array' works for undirected networks", {
  mplearray <- ergmMPLE(flomarriage~edges+absdiff("wealth"), output="array")
  ut <- upper.tri(mplearray$response)
  lt <- lower.tri(mplearray$response,diag=TRUE)
  nas <- rep(NA_real_,network.dyadcount(flomarriage)+network.size(flomarriage))
  ones <- rep(1,network.dyadcount(flomarriage))
  zeros <- numeric(network.dyadcount(flomarriage)+network.size(flomarriage))

  expect_equal(mplearray$response[ut], as.matrix(flomarriage)[ut], ignore_attr=TRUE)
  expect_equal(mplearray$response[lt], nas, ignore_attr=TRUE)

  expect_equal(mplearray$predictor[,,1][ut], ones, ignore_attr=TRUE)
  expect_equal(mplearray$predictor[,,1][lt], nas, ignore_attr=TRUE)

  wealth <- flomarriage%v%"wealth"
  expect_equal(mplearray$predictor[,,2][ut], outer(wealth, wealth, FUN=function(x,y) abs(x-y))[ut], ignore_attr=TRUE)
  expect_equal(mplearray$predictor[,,2][lt], nas, ignore_attr=TRUE)

  expect_equal(mplearray$weights[ut], ones, ignore_attr=TRUE)
  expect_equal(mplearray$weights[lt], zeros, ignore_attr=TRUE)
})

bfl <- get.inducedSubgraph(flomarriage, 1:7, 8:16)

test_that("output = 'array' works for bipartite networks with expand.bipartite = FALSE", {
  mplearray <- ergmMPLE(bfl~edges+absdiff("wealth"), output="array")
  ones <- rep(1,network.dyadcount(bfl))

  expect_equal(mplearray$response, as.matrix(bfl), ignore_attr=TRUE)

  expect_equal(mplearray$predictor[,,1], ones, ignore_attr=TRUE)

  wealth <- flomarriage%v%"wealth"
  expect_equal(mplearray$predictor[,,2], outer(wealth[1:7], wealth[8:16], FUN=function(x,y) abs(x-y)), ignore_attr=TRUE)

  expect_equal(mplearray$weights, ones, ignore_attr=TRUE)
})


test_that("output = 'array' works for bipartite networks with expand.bipartite = TRUE", {
  mplearray <- ergmMPLE(bfl~edges+absdiff("wealth"), output="array", expand.bipartite=TRUE)
  ut <- upper.tri(mplearray$response)
  lt <- lower.tri(mplearray$response,diag=TRUE)
  b <- rep(c(TRUE,FALSE), c(7,9))
  b <- outer(b,b, `!=`)
  ut <- ut&b
  lt <- lt|!b

  nas <- rep(NA_real_,sum(lt))
  ones <- rep(1,sum(ut))
  zeros <- numeric(sum(lt))

  expect_equal(mplearray$response[ut], as.matrix(flomarriage)[ut], ignore_attr=TRUE)
  expect_equal(mplearray$response[lt], nas, ignore_attr=TRUE)

  expect_equal(mplearray$predictor[,,1][ut], ones, ignore_attr=TRUE)
  expect_equal(mplearray$predictor[,,1][lt], nas, ignore_attr=TRUE)

  wealth <- flomarriage%v%"wealth"
  expect_equal(mplearray$predictor[,,2][ut], outer(wealth, wealth, FUN=function(x,y) abs(x-y))[ut], ignore_attr=TRUE)
  expect_equal(mplearray$predictor[,,2][lt], nas, ignore_attr=TRUE)

  expect_equal(mplearray$weights[ut], ones, ignore_attr=TRUE)
  expect_equal(mplearray$weights[lt], zeros, ignore_attr=TRUE)
})
