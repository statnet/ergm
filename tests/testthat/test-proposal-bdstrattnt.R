#  File tests/testthat/test-proposal-bdstrattnt.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################


test_that("BDStratTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(1000, 125, 125, 350)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  
  
  ## redo above with lower target stats and a nontrivial upper bound on degree
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(160, 20, 20, 100)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(530, 30, 30, 450)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))  
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  
  ## again but now also with an levels2
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))
  nw %v% "sex" <- rep(c("M","F"), 500)
  
  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(80, 10, 10, 50)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
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
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
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
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(100, 10, 0, 0, 100, 100, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd + Strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo above but with a nontrivial upper bound on degree
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 90, 45, 0, 45, 2, 2, 0)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + Strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  
  
  # ditto but also with levels2
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), 0)
  
  # this time with only same-sex ties
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + Strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), network.edgecount(nws2))
})

test_that("BDStratTNT works with churning", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 30), rep("B", 30), rep("W", 940))
  nw %v% "sex"  <- rep(c("W", "X", "Y", "Z"), 250)
  
  levels2 <- matrix(0, 4, 4)
  levels2[1,3] <- levels2[3,1] <- levels2[2,2] <- levels2[3,4] <- levels2[4,3] <- levels2[4,4] <- 1
  levels2 <- levels2 > 0
  
  pmat <- matrix(c(25, 50, 5, 50, 25, 5, 5, 5, 100),3,3,byrow=TRUE)

  # impossible to hit these exactly
  target.stats <- c(211, 25, 50, 25, 5, 5, 100)
  nws <- san(nw ~ edges + nodemix("race",levels2=TRUE), target.stats = target.stats, constraints = ~bd(maxout=4) + blocks(attr="sex", levels2=levels2) + Strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats + 1))
  expect_equal(unname(summary(nws ~ degrange(5))), 0)
  # and check sex nodemix
  srs <- summary(nws ~ nodemix("sex",levels2=TRUE))
  expect_true(all(srs[as.logical(levels2[upper.tri(levels2,diag=TRUE)])] == 0))
})

test_that("BDStratTNT simulates reasonably", {
  for(deg_bound in 1:5) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c("X","Y","Z"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    
    levels2 <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
          
    nw_sim <- nw
    
    for(i in 1:5) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "sex", levels2 = levels2) + Strat(attr = "vattr", pmat = pmat),
                         output = "network")
      summ_stats <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE) + nodemix("sex",levels2=TRUE) + degrange(deg_bound + 1))
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
})

test_that("BDStratTNT works with degree bound saturation", {
  nw <- network.initialize(900, dir=FALSE)

  nw %v% "race" <- rep(c("A","B","C"), times = c(30,30,840))
  nw %v% "sex" <- rep(c("X","Y","Z"), length.out=900)

  pmat <- matrix(1,3,3)

  # mix.race.A.A mix.race.A.B mix.race.B.B mix.race.A.C mix.race.B.C mix.race.C.C

  target.stats <- c(425, 10, 5, 10, 0, 0, 400.01)
  nws <- san(nw ~ edges + nodemix("race",levels2=TRUE), target.stats = target.stats, constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = matrix(c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE),3,3)) + Strat(attr = "race", pmat = pmat), control=control.san(SAN.invcov.diag=TRUE, SAN.maxit = 4, SAN.nsteps=5e4))
  sr <- summary(nws ~ edges + nodemix("race",levels2=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= pmax(1, 0.05*target.stats)))
  expect_true(all(summary(nws ~ concurrent + nodemix("sex", levels2=c(1,5))) == 0))
})

test_that("BDStratTNT constrains undirected appropriately", {
  nw <- network.initialize(100, dir=FALSE)
  nw %v% "attr" <- rep(c("A","B","C","D","E"), each = 20)
  nw %v% "strat_attr" <- rep(1:3, length.out=100)
  nw[cbind(1:10,30:21)] <- 1
  nw[cbind(44:53,99:90)] <- 1
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~attr, levels2=c(2,13)) + Strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,30:21)] == 1))
  expect_true(all(nws[cbind(44:53,99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2=13)) == 10)
  expect_true(summary(nws ~ edges) > 1000)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1) + blocks(~attr, levels2=c(2,13)) + Strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,30:21)] == 1))
  expect_true(all(nws[cbind(44:53,99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2=13)) == 10)
  expect_true(summary(nws ~ edges) > 30)  
})

test_that("BDStratTNT constrains bipartite appropriately", {
  nw <- network.initialize(100, bip=30, dir=FALSE)
  nw %v% "attr" <- c(rep(c("A","B","C"), each = 10), rep(c("D","E","F","G"), times = c(20,20,20,10)))
  nw %v% "strat_attr" <- rep(1:3, length.out=100)
  nw[cbind(1:10,100:91)] <- 1
  nw[cbind(25:21,60:56)] <- 1
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~attr, levels2=c(6,10)) + Strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,100:91)] == 1))
  expect_true(all(nws[cbind(25:21,60:56)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=6)) == 5)
  expect_true(summary(nws ~ nodemix(~attr, levels2=10)) == 10)
  expect_true(summary(nws ~ edges) > 500)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1) + blocks(~attr, levels2=c(6,10)) + Strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,100:91)] == 1))
  expect_true(all(nws[cbind(25:21,60:56)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=6)) == 5)
  expect_true(summary(nws ~ nodemix(~attr, levels2=10)) == 10)
  expect_true(summary(nws ~ edges) > 20)
})

test_that("BDStratTNT handles undirected arguments correctly", {
  nw <- network.initialize(100, dir=FALSE)
  nw %v% "bd_attr" <- rep(1:3, length.out=100)
  nw %v% "strat_attr" <- rep(1:7, length.out=100)
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1) + Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1) + Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
})

test_that("BDStratTNT handles bipartite arguments correctly", {
  nw <- network.initialize(100, dir=FALSE, bip=30)
  nw %v% "bd_attr" <- c(rep(1:3, length.out=30), rep(6:10, length.out=70))
  nw %v% "strat_attr" <- rep(1:7, length.out=100)

  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + bd(maxout = 1), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~Strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + bd(maxout = 1) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
})

