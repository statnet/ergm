#  File tests/testthat/test-proposal-strattnt.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################


test_that("StratTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=1e4))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(1000, 125, 125, 350)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=2e4))
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})


test_that("StratTNT works with directed networks", {
  nw <- network.initialize(1000, dir=TRUE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(c(100, 350, 0, 10, 100, 0, 100, 0, 840),3,3,byrow=TRUE)

  target.stats <- c(100, 10, 100, 350, 100, 0, 0, 0, 840)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=1e4))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(50, 50, 350, 50, 50, 100, 50, 400, 400),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat3, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=2e4))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})

test_that("StratTNT works with bipartite networks", {
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 100, 2, 2, 0, 2, 100, 100, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 100, 100, 0, 100, 2, 2, 0)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=1e4))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(100, 10, 0, 0, 100, 100, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat3, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=2e4))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})

test_that("StratTNT works with churning", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 30), rep("B", 30), rep("W", 940))

  pmat <- matrix(c(50, 50, 5, 50, 50, 5, 5, 5, 100),3,3,byrow=TRUE)

  # impossible to hit these exactly
  target.stats <- c(261, 50, 50, 50, 5, 5, 100)
  nws <- san(nw ~ edges + nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~strat(pmat=pmat, attr="race"))
  sr <- summary(nws ~ edges + nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats + 1))
})

test_that("StratTNT simulates reasonably", {

  net_size <- 500L

  nw <- network.initialize(net_size, dir = FALSE)

  vattr <- sample(c("A","B","C"), net_size, TRUE)
  
  nw %v% "vattr" <- vattr
  
  pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
      
  nw_sim <- nw
  
  for(i in 1:5) {
    nw_sim <- simulate(nw_sim ~ edges, 
                       coef = c(-3), 
                       constraints = ~strat(attr = "vattr", pmat = pmat),
                       output = "network")
    summ_stats <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
    expect_true(summ_stats["mix.vattr.A.A"] == 0)
    expect_true(summ_stats["mix.vattr.B.B"] == 0)
    expect_true(summ_stats["mix.vattr.A.B"] > 0)
    expect_true(summ_stats["mix.vattr.A.C"] > 0)
    expect_true(summ_stats["mix.vattr.B.C"] > 0)
    expect_true(summ_stats["mix.vattr.C.C"] > 0)    
  }  
})

test_that("BDStratTNT handles undirected arguments correctly", {
  nw <- network.initialize(100, dir=FALSE)
  nw %v% "strat_attr" <- rep(1:3, length.out=100)

  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(0,rep(1,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-1)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(~strat_attr, pmat = matrix(c(0,rep(1,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(1,0,1,0,rep(1,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-2)) > 0))
})

test_that("BDStratTNT handles directed arguments correctly", {
  nw <- network.initialize(100, dir=TRUE)
  nw %v% "strat_attr" <- rep(1:3, length.out=100)

  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(0,rep(1,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-1)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(~strat_attr, pmat = matrix(c(0,rep(1,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(1,0,1,1,rep(1,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-2)) > 0))
})

test_that("BDStratTNT handles bipartite arguments correctly", {
  nw <- network.initialize(100, dir=FALSE, bip=30)
  nw %v% "strat_attr" <- c(rep(1:3, length.out=30), rep(6:10, length.out=70))

  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(0,rep(1,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-1)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(~strat_attr, pmat = matrix(c(0,rep(1,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(strat_attr = ~strat_attr, pmat = matrix(c(1,0,1,1,rep(1,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~strat_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~strat_attr, levels2=-2)) > 0))
})

