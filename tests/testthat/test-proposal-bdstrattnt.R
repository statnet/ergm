#  File tests/testthat/test-proposal-bdstrattnt.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

context("test-proposal-bdstrattnt.R")

test_that("BDStratTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(1000, 125, 125, 350)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  
  
  ## redo above with lower target stats and a nontrivial upper bound on degree
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(160, 20, 20, 100)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(530, 30, 30, 450)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))  
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  
  ## again but now also with an fmat
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))
  nw %v% "sex" <- rep(c("M","F"), 500)
  
  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(80, 10, 10, 50)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))  
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), 0)

  ## this time with only same-sex ties
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))
  nw %v% "sex" <- rep(c("M","F"), 500)
  
  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(80, 10, 10, 50)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))  
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), network.edgecount(nws2))

  
})

test_that("BDStratTNT works with bipartite networks", {
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 100, 2, 2, 0, 2, 100, 100, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 100, 100, 0, 100, 2, 2, 0)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(100, 10, 0, 0, 100, 100, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(pmat=pmat3, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo above but with a nontrivial upper bound on degree
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 90, 45, 0, 45, 2, 2, 0)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat3, Strat_attr="race")), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  
  
  # ditto but also with fmat
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(1,0,0,1), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat3, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(1,0,0,1), 2 ,2))), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), 0)
  
  # this time with only same-sex ties
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3, SAN.prop.args = list(bound = 5, pmat=pmat, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(0,1,1,0), 2, 2))), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race"), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4, SAN.prop.args = list(bound = 5, pmat=pmat3, Strat_attr="race", BD_attr = "sex", fmat = matrix(c(0,1,1,0), 2 ,2))), constraints="BDStratTNT"~.)
  sr <- summary(nws2 ~ nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), network.edgecount(nws2))
})

test_that("BDStratTNT works with churning", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 30), rep("B", 30), rep("W", 940))
  nw %v% "sex"  <- rep(c("W", "X", "Y", "Z"), 250)
  
  fmat <- matrix(0, 4, 4)
  fmat[1,3] <- fmat[3,1] <- fmat[2,2] <- fmat[3,4] <- fmat[4,3] <- fmat[4,4] <- 1
  
  pmat <- matrix(c(25, 50, 5, 50, 25, 5, 5, 5, 100),3,3,byrow=TRUE)

  # impossible to hit these exactly
  target.stats <- c(211, 25, 50, 25, 5, 5, 100)
  nws <- san(nw ~ edges + nodemix("race"), target.stats = target.stats, control=control.san(SAN.prop.args = list(bound = 4, BD_attr = "sex", fmat = fmat, Strat_attr = "race", pmat = pmat)), constraints="BDStratTNT"~.)
  sr <- summary(nws ~ edges + nodemix("race"))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats + 1))
  expect_equal(unname(summary(nws ~ degrange(5))), 0)
  # and check sex nodemix
  srs <- summary(nws ~ nodemix("sex"))
  expect_true(all(srs[as.logical(fmat[upper.tri(fmat,diag=TRUE)])] == 0))
})

test_that("BDStratTNT simulates reasonably", {
  for(deg_bound in 1:5) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c("X","Y","Z"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    
    fmat <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
      
    control <- control.simulate.formula(MCMC.prop.weights = "BDStratTNT", 
                                        MCMC.prop.args = list(bound = deg_bound, 
                                                              BD_attr = "sex",
                                                              fmat = fmat,
                                                              Strat_attr = "vattr",
                                                              pmat = pmat))
    
    nw_sim <- nw
    
    for(i in 1:5) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         output = "network",
                         control = control)
      summ_stats <- summary(nw_sim ~ nodemix("vattr") + nodemix("sex") + degrange(deg_bound + 1))
      expect_true(summ_stats["mix.vattr.A.A"] == 0)
      expect_true(summ_stats["mix.vattr.B.B"] == 0)
      expect_true(summ_stats["mix.vattr.A.B"] > 0)
      expect_true(summ_stats["mix.vattr.A.C"] > 0)
      expect_true(summ_stats["mix.vattr.B.C"] > 0)
      expect_true(summ_stats["mix.vattr.C.C"] > 0)    
  
      expect_true(summ_stats["mix.sex.X.X"] == 0)
      expect_true(summ_stats["mix.sex.X.Z"] == 0)
      expect_true(summ_stats["mix.sex.X.Y"] > 0)
      expect_true(summ_stats["mix.sex.Y.Y"] > 0)
      expect_true(summ_stats["mix.sex.Y.Z"] > 0)
      expect_true(summ_stats["mix.sex.Z.Z"] > 0)
  
      expect_true(summ_stats[paste0("deg", deg_bound + 1, "+")] == 0)    
    }  
  }
}
