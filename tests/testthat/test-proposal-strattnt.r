#  File tests/testthat/test-proposal-strattnt.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

context("test-proposal-strattnt.R")

test_that("StratTNT works with undirected unipartite networks", {
  nw <- network(1000, dir=FALSE, numedges=0)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.nsteps=4e4, SAN.prop.args = list(pmat=pmat, attr="race")), constraints="StratTNT"~.)
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})


test_that("StratTNT works with directed networks", {
  nw <- network(1000, dir=TRUE, numedges=0)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(c(100, 350, 0, 10, 100, 0, 100, 0, 840),3,3,byrow=TRUE)

  target.stats <- c(100, 10, 100, 350, 100, 0, 0, 0, 840)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.nsteps=5e5, SAN.prop.args = list(pmat=pmat, attr="race")), constraints="StratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})

test_that("StratTNT works with bipartite networks", {
  nw <- network(900, bip = 100, numedges=0)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 100, 2, 2, 0, 2, 100, 100, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 100, 100, 0, 100, 2, 2, 0)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.nsteps=5e5, SAN.prop.args = list(pmat=pmat, attr="race")), constraints="StratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))  
})
