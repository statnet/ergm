#  File tests/testthat/test-proposal-bdtnt.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

context("test-proposal-bdtnt.R")

test_that("BDTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  target.stats <- c(500)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=2e3, SAN.prop.args = list(bound = 1)), constraints="BDTNT"~.)
  sr <- summary(nws ~ edges + concurrent)  
  
  expect_equal(unname(sr), c(500,0))
  
  
  
  target.stats <- c(1000)
  nws2 <- san(nws ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=2e3, SAN.prop.args = list(bound = 2)), constraints="BDTNT"~.)
  sr2 <- summary(nws2 ~ edges + degree(2) + degrange(3))  
  
  expect_true(all(abs(sr2 - c(1000,1000, 0)) <= c(1,2,0)))
  
  
  
  target.stats <- c(1500)
  nws22 <- san(nws2 ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=2e3, SAN.prop.args = list(bound = 3)), constraints="BDTNT"~.)
  sr22 <- summary(nws22 ~ edges + degree(3) + degrange(4))  
  
  expect_true(all(abs(sr22 - c(1500,1000, 0)) <= c(2,2,0)))
  
  

  ## may be off by small amount  
  target.stats <- c(1000)
  nws2a <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e3, SAN.prop.args = list(bound = 2)), constraints="BDTNT"~.)
  sr2a <- summary(nws2a ~ edges + degree(2) + degrange(3))  
  
  expect_true(all(abs(sr2a - c(1000,1000, 0)) <= c(1,2,0)))
  
  
  ## may be off by small amount  
  target.stats <- c(1500)
  nws22a <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=6e3, SAN.prop.args = list(bound = 3)), constraints="BDTNT"~.)
  sr22a <- summary(nws22a ~ edges + degree(3) + degrange(4))  
  
  expect_true(all(abs(sr22a - c(1500,1000, 0)) <= c(2,2,0)))
  
    
  
  
  nw %v% "sex" <- rep(c("A","B"), 500)
  target.stats <- c(500)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=2e3, SAN.prop.args = list(bound = 1, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr3 <- summary(nws ~ edges + concurrent + nodematch("sex"))  
  
  expect_equal(unname(sr3), c(500,0, 0))
  

  target.stats <- c(500)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=2e3, SAN.prop.args = list(bound = 1, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr4 <- summary(nws ~ edges + concurrent + nodematch("sex"))  
  
  expect_equal(unname(sr4), c(500,0, 500))
  
  ## may be off by small amount
  target.stats <- c(1500)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=6e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr5 <- summary(nws ~ edges + degree(3) + degrange(4) + nodematch("sex"))  
  
  expect_true(all(abs(sr5 - c(1500, 1000, 0, 0)) <= c(2,4,0,2)))
  
  ## may be off by small amount
  target.stats <- c(1500)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=6e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr6 <- summary(nws ~ edges + degree(3) + degrange(4) + nodematch("sex"))  
  
  expect_true(all(abs(sr6 - c(1500, 1000, 0, 1500)) <= c(4,4,0,4)))
  

  target.stats <- c(1000)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr7 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))  
  
  expect_equal(unname(sr7), c(1000,0, 0))
  

  target.stats <- c(1000)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr8 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))  
  
  expect_equal(unname(sr8), c(1000,0, 1000))


})

test_that("BDTNT works with bipartite networks", {
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  target.stats <- c(100)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e2, SAN.prop.args = list(bound = 1)), constraints="BDTNT"~.)
  sr <- summary(nws ~ edges + b1degree(1) + degrange(2))  
  
  expect_equal(unname(sr), c(100,100, 0))
  
  target.stats <- c(200)
  nws2 <- san(nws ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e2, SAN.prop.args = list(bound = 2)), constraints="BDTNT"~.)
  sr2 <- summary(nws2 ~ edges + b1degree(2) + degrange(3))  
  
  expect_true(all(abs(sr2 - c(200,100, 0)) <= c(0,0,0)))
  
  target.stats <- c(300)
  nws22 <- san(nws2 ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e2, SAN.prop.args = list(bound = 3)), constraints="BDTNT"~.)
  sr22 <- summary(nws22 ~ edges + b1degree(3) + degrange(4))  
  
  expect_true(all(abs(sr22 - c(300,100, 0)) <= c(0,0,0)))
  
  target.stats <- c(200)
  nws2a <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=8e2, SAN.prop.args = list(bound = 2)), constraints="BDTNT"~.)
  sr2a <- summary(nws2a ~ edges + b1degree(2) + degrange(3))  
  
  expect_true(all(abs(sr2a - c(200,100, 0)) <= c(0,0,0)))
  
  target.stats <- c(300)
  nws22a <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1.2e3, SAN.prop.args = list(bound = 3)), constraints="BDTNT"~.)
  sr22a <- summary(nws22a ~ edges + b1degree(3) + degrange(4))  
  
  expect_true(all(abs(sr22a - c(300,100, 0)) <= c(0,0,0)))
  
    
  
  
  nw %v% "sex" <- c(rep(c("A","B"), 50), rep(c("A","B"), 450))
  target.stats <- c(100)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e2, SAN.prop.args = list(bound = 1, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr3 <- summary(nws ~ edges + concurrent + nodematch("sex"))  
  
  expect_equal(unname(sr3), c(100,0, 0))
  
  target.stats <- c(100)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=4e2, SAN.prop.args = list(bound = 1, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr4 <- summary(nws ~ edges + concurrent + nodematch("sex"))  
  
  expect_equal(unname(sr4), c(100,0, 100))
  
  target.stats <- c(300)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1.2e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr5 <- summary(nws ~ edges + b1degree(3) + degrange(4) + nodematch("sex"))  
  
  expect_true(all(abs(sr5 - c(300, 100, 0, 0)) <= c(0,0,0,0)))
  
  target.stats <- c(300)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1.2e3, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr6 <- summary(nws ~ edges + b1degree(3) + degrange(4) + nodematch("sex"))  
  
  expect_true(all(abs(sr6 - c(300, 100, 0, 300)) <= c(0,0,0,0)))
  

  target.stats <- c(200)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=8e2, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDTNT"~.)
  sr7 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))  
  
  expect_equal(unname(sr7), c(200,0, 0))
  

  target.stats <- c(200)
  nws <- san(nw ~ edges , target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=8e2, SAN.prop.args = list(bound = 3, attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDTNT"~.)
  sr8 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))  
  
  expect_equal(unname(sr8), c(200,0, 200))

})
