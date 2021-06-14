#  File tests/testthat/test-proposal-bdstrattnt.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################


test_that("BDStratTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(1000, 125, 125, 350)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd + strat(attr = "race", pmat = pmat))
  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  
  
  ## redo above with lower target stats and a nontrivial upper bound on degree
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(1,3,3)
  diag(pmat) <- c(2,2,10)

  target.stats <- c(160, 20, 20, 100)
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(530, 30, 30, 450)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + strat(attr = "race", pmat = pmat))
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
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
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
  nws <- san(nw ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE))  
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10,3,3)
  diag(pmat) <- c(7,7,20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race",levels=NULL, diff=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
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
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(100, 10, 0, 0, 100, 100, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd + strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # redo above but with a nontrivial upper bound on degree
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- c(0, 2, 90, 45, 0, 45, 2, 2, 0)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + strat(attr = "race", pmat = pmat3))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  
  
  # ditto but also with levels2
  nw <- network.initialize(900, bip = 100, dir=FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W",860))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0),3,3,byrow=TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=diag(TRUE, 2)) + strat(attr = "race", pmat = pmat3))
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
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=5e3), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + strat(attr = "race", pmat = pmat))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))
  
  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, control=control.san(SAN.maxit = 1, SAN.nsteps=1e4), constraints = ~bd(maxout=5) + blocks(attr="sex", levels2=!diag(TRUE, 2)) + strat(attr = "race", pmat = pmat3))
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
  nws <- san(nw ~ edges + nodemix("race",levels2=TRUE), target.stats = target.stats, constraints = ~bd(maxout=4) + blocks(attr="sex", levels2=levels2) + strat(attr = "race", pmat = pmat))
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
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "vattr", pmat = pmat),
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

test_that("BDStratTNT simulates reasonably with heterogeneous degree bounds", {
  for(deg_bound in 1:5) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c(1,2,3), net_size, TRUE)
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)
    
    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(1,2,2,4), c(5,2,3,4))] <- TRUE
    blocks_levels_2 <- blocks_levels_2 | t(blocks_levels_2)
    
    levels2 <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
          
    nw_sim <- nw

    maxout <- matrix(0, nrow = net_size, ncol = 3)

    for(row_index in 1:3) {
      for(col_index in 1:3) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    for(i in 1:5) {    
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         constraints = ~bd(attribs = attribs, maxout = maxout) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "network")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,3)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,3)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3,5,10,11)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3,5,10,11)] > 0))
      
      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)), to = factor(sex[c(el[,c(2,1)])], levels = seq_len(3)))
      expect_true(all(degs <= maxout))
    }  
  }
})

skip("Skipping the rest for time.")

test_that("BDStratTNT simulates reasonably with bipartite heterogeneous degree bounds", {
  for(deg_bound in 1:5) {
    net_size <- 2000L
    bip <- 700L
    
    nw <- network.initialize(net_size, dir = FALSE, bip = bip)
  
    vattr <- c(sample(c("A","B","C","D"), bip, TRUE), sample(c("X","Y","Z"), net_size - bip, TRUE))
    sex <- c(sample(c(1,2,3,4,5), bip, TRUE), sample(c(6,7,8,9,10,11), net_size - bip, TRUE))
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = length(unique(sex)))
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- c(sample(c(1,2,3), bip, TRUE), sample(c(4,5,6,7), net_size - bip, TRUE))
    
    blocks_levels_2 <- matrix(FALSE, nrow = 3, 4)
    blocks_levels_2[cbind(c(3,2,2), c(1,2,3))] <- TRUE
    
    levels2 <- matrix(as.logical(round(runif(11*11))), nrow = 11, ncol = 11)
    levels2 <- levels2 | t(levels2)
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),nrow = 4, ncol = 3)
          
    nw_sim <- nw

    maxout <- matrix(0, nrow = net_size, ncol = 11)

    for(row_index in 1:11) {
      for(col_index in 1:11) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    for(i in 1:5) {    
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         constraints = ~bd(attribs = attribs, maxout = maxout) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "network")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,5,7,12)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,5,7,12)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3,5,8)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3,5,8)] > 0))
      
      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)), to = factor(sex[c(el[,c(2,1)])], levels = seq_len(11)))
      expect_true(all(degs <= maxout))
    }  
  }
})

test_that("BDStratTNT simulates reasonably with directed heterogeneous degree bounds", {
  for(deg_bound in 1:5) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = TRUE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c(1,2,3), net_size, TRUE)
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)
    
    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(5,2,2,4), c(1,2,3,4))] <- TRUE
    
    levels2 <- matrix(c(1,0,0,0,0,1,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,0,0,0,1,0),3,3)
          
    nw_sim <- nw
    
    maxout <- matrix(0, nrow = net_size, ncol = 3)

    for(row_index in 1:3) {
      for(col_index in 1:3) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    maxin <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxin[maxin < 0] <- 0
    
    for(i in 1:5) {      
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         constraints = ~bd(attribs = attribs, maxout = maxout, maxin = maxin) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "network")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,8)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,8)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(5,8,14,22)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(5,8,14,22)] > 0))
      
      el <- as.edgelist(nw_sim)
      out_degs <- table(from = factor(c(el[,1]), levels = seq_len(net_size)), to = factor(sex[c(el[,2])], levels = seq_len(3)))
      in_degs <- table(from = factor(c(el[,2]), levels = seq_len(net_size)), to = factor(sex[c(el[,1])], levels = seq_len(3)))
      expect_true(all(out_degs <= maxout))
      expect_true(all(in_degs <= maxin))      
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
  nws <- san(nw ~ edges + nodemix("race",levels2=TRUE), target.stats = target.stats, constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = matrix(c(TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE),3,3)) + strat(attr = "race", pmat = pmat), control=control.san(SAN.invcov.diag=TRUE, SAN.maxit = 4, SAN.nsteps=5e4))
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
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~attr, levels2=c(2,13)) + strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,30:21)] == 1))
  expect_true(all(nws[cbind(44:53,99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2=13)) == 10)
  expect_true(summary(nws ~ edges) > 1000)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1) + blocks(~attr, levels2=c(2,13)) + strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,30:21)] == 1))
  expect_true(all(nws[cbind(44:53,99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2=13)) == 10)
  expect_true(summary(nws ~ edges) > 30)  

  nw <- network.initialize(100, dir=FALSE)
  nw %v% "attr" <- rep(c("B","A","C","D","E"), each = 20)
  nw %v% "strat_attr" <- rep(1:3, length.out=100)
  nw[cbind(1:10,30:21)] <- 1
  nw[cbind(44:53,99:90)] <- 1
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~attr, levels2=c(2,13)) + strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,30:21)] == 1))
  expect_true(all(nws[cbind(44:53,99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2=13)) == 10)
  expect_true(summary(nws ~ edges) > 1000)
})

test_that("BDStratTNT constrains bipartite appropriately", {
  nw <- network.initialize(100, bip=30, dir=FALSE)
  nw %v% "attr" <- c(rep(c("A","B","C"), each = 10), rep(c("D","E","F","G"), times = c(20,20,20,10)))
  nw %v% "strat_attr" <- rep(1:3, length.out=100)
  nw[cbind(1:10,100:91)] <- 1
  nw[cbind(25:21,60:56)] <- 1
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~attr, levels2=c(6,10)) + strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
  expect_true(all(nws[cbind(1:10,100:91)] == 1))
  expect_true(all(nws[cbind(25:21,60:56)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2=6)) == 5)
  expect_true(summary(nws ~ nodemix(~attr, levels2=10)) == 10)
  expect_true(summary(nws ~ edges) > 500)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout = 1) + blocks(~attr, levels2=c(6,10)) + strat(~strat_attr, pmat = matrix(2 + runif(9),3,3)))
  
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

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1) + strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~bd(maxout=1) + strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,8)),3,3)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,TRUE,rep(FALSE,5)),3,3))))
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

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-1)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + bd(maxout = 1), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
  
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))

  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)) + bd(maxout = 1) + blocks(~bd_attr, levels2 = matrix(c(TRUE,rep(FALSE,14)),nrow=3,ncol=5)), control = list(MCMC.prop.weights = "BDStratTNT", MCMC.prop.args = list(blocks_attr = ~bd_attr, levels2 = matrix(c(FALSE,TRUE,FALSE,FALSE,rep(FALSE,11)),nrow=3,ncol=5))))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)  
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels=TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
})

test_that("BDStratTNT handles atypical levels specifications correctly", {
  nw <- network.initialize(100, dir=FALSE)
  nw %v% "bd_attr" <- rep(1:3, length.out=100)
  nw %v% "strat_attr" <- rep(1:5, length.out = 100)
  pmat <- matrix(2 + runif(25), 5, 5)
  
  ## should be unconstrained
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))

  ## should also be unconstrained
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels=FALSE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))

  ## any pairing with a 3 should be allowed, with all other pairings forbidden
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels=I(c(1,2,4,6)), levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=c(4,5,6))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-c(4,5,6))) == 0))

  ## only 2-2 pairings should be allowed
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels=I(c(1,2,3,4,6)), levels2=-3) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=c(3))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-c(3))) == 0))
  
  ## should fail as we omit all pairings
  expect_error(nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat)))
  
  
  ## similar bipartite tests
  nw <- network.initialize(100, dir=FALSE, bip = 30)
  nw %v% "bd_attr" <- c(rep(1:3, length.out=30), rep(10:16, length.out = 70))
  nw %v% "strat_attr" <- c(rep(1:5, length.out = 30), rep(1:4, length.out = 70))
  pmat <- matrix(2 + runif(20), nrow = 5, ncol = 4)

  ## should be unconstrained
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=TRUE, b2levels=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))

  ## should also be unconstrained
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=FALSE, b2levels=FALSE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=FALSE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b2levels=FALSE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=FALSE, b2levels=FALSE, levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=FALSE, levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b2levels=FALSE, levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=TRUE)) > 0))

  ## any pairing with a 3 should be allowed, with all other pairings forbidden
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=I(c(1,2,4,6)), levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=3*(1:7))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-3*(1:7))) == 0))

  ## only 1-14 pairings should be allowed
  nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, b1levels=I(c(1,2,3,4,6)), levels2=-21) + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, b1levels = I(1), b2levels = I(14), levels2=TRUE)) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-c(13))) == 0))
  
  ## should fail as we omit all pairings
  expect_error(nws <- simulate(nw ~ edges, coef = c(0), constraints = ~blocks(~bd_attr, levels2=TRUE) + strat(attr = ~strat_attr, pmat = pmat)))
})


test_that("BDStratTNT works with directed networks", {
  nw <- network.initialize(1000, dir=TRUE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W",960))

  pmat <- matrix(c(100, 350, 0, 10, 100, 0, 100, 0, 840),3,3,byrow=TRUE)

  target.stats <- c(100, 10, 100, 350, 100, 0, 0, 0, 840)
  nws <- san(nw ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~bd(maxout=40, maxin=40) + strat(pmat=pmat, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=1e4))
  sr <- summary(nws ~ nodemix("race",levels2=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(50, 50, 350, 50, 50, 100, 50, 400, 400),3,3,byrow=TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)
  nws2 <- san(nws ~ nodemix("race",levels2=TRUE), target.stats = target.stats, constraints=~bd(maxout=40, maxin=40) + strat(pmat=pmat3, attr="race"), control=control.san(SAN.maxit = 1, SAN.nsteps=2e4))
  sr <- summary(nws2 ~ nodemix("race",levels2=TRUE))
  
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})

test_that("BDStratTNT simulates directed reasonably", {

  net_size <- 1000L

  nw <- network.initialize(net_size, dir = TRUE)

  vattr <- sample(c("A","B","C"), net_size, TRUE)
  
  nw %v% "vattr" <- vattr
  nw %v% "sex" <- sample(c("X","Y","Z"), net_size, TRUE)
  
  pmat <- 1 - matrix(c(1,0,0,1,1,0,0,1,0),3,3)
      
  nw_sim <- nw
  
  for(i in 1:5) {
    nw_sim <- simulate(nw_sim ~ edges, 
                       coef = c(0), 
                       constraints = ~bd(maxout = i, maxin = i + 1) + blocks(attr = "sex", levels2 = matrix(c(TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE),3,3)) + strat(attr = "vattr", pmat = pmat),
                       output = "network")
    summ_stats <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE) + nodemix("sex", levels2 = TRUE) + odegrange(i + 1) + idegrange(i + 2))
    expect_true(summ_stats["mix.vattr.A.A"] == 0)
    expect_true(summ_stats["mix.vattr.B.B"] == 0)
    expect_true(summ_stats["mix.vattr.A.B"] == 0)
    expect_true(summ_stats["mix.vattr.B.C"] == 0)
    expect_true(summ_stats["mix.vattr.A.C"] > 0)
    expect_true(summ_stats["mix.vattr.B.A"] > 0)
    expect_true(summ_stats["mix.vattr.C.A"] > 0)
    expect_true(summ_stats["mix.vattr.C.B"] > 0)    
    expect_true(summ_stats["mix.vattr.C.C"] > 0)    
    expect_true(summ_stats["mix.sex.X.X"] == 0)
    expect_true(summ_stats["mix.sex.X.Y"] > 0)
    expect_true(summ_stats["mix.sex.X.Z"] > 0)
    expect_true(summ_stats["mix.sex.Y.X"] > 0)
    expect_true(summ_stats["mix.sex.Y.Y"] > 0)
    expect_true(summ_stats["mix.sex.Y.Z"] == 0)
    expect_true(summ_stats["mix.sex.Z.X"] == 0)
    expect_true(summ_stats["mix.sex.Z.Y"] == 0)    
    expect_true(summ_stats["mix.sex.Z.Z"] > 0)    
    expect_true(summ_stats[paste0("odeg", i + 1, "+")] == 0)
    expect_true(summ_stats[paste0("ideg", i + 2, "+")] == 0)    
  }  
})
