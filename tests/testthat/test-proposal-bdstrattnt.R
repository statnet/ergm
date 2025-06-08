#  File tests/testthat/test-proposal-bdstrattnt.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

test_that("BDStratTNT works with undirected unipartite networks and bd/blocks only", {
  nw <- network.initialize(1000, directed = FALSE)

  target.stats <- c(500)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 2e3))

  sr <- summary(nws ~ edges + concurrent)
  expect_equal(unname(sr), c(500, 0))

  target.stats <- c(1000)
  nws2 <- san(nws ~ edges,
              target.stats = target.stats,
              constraints = ~bd(maxout = 2),
              control = control.san(SAN.maxit = 1, SAN.nsteps = 2e3))

  sr2 <- summary(nws2 ~ edges + degree(2) + degrange(3))
  expect_true(all(abs(sr2 - c(1000, 1000, 0)) <= c(1, 2, 0)))

  target.stats <- c(1500)
  nws22 <- san(nws2 ~ edges,
               target.stats = target.stats,
               constraints = ~bd(maxout=3),
               control = control.san(SAN.maxit = 1, SAN.nsteps = 2e3))

  sr22 <- summary(nws22 ~ edges + degree(3) + degrange(4))
  expect_true(all(abs(sr22 - c(1500, 1000, 0)) <= c(2, 2, 0)))

  ## may be off by small amount
  target.stats <- c(1000)
  nws2a <- san(nw ~ edges,
               target.stats = target.stats,
               constraints = ~bd(maxout = 2),
               control=control.san(SAN.maxit = 1, SAN.nsteps = 4e3))

  sr2a <- summary(nws2a ~ edges + degree(2) + degrange(3))
  expect_true(all(abs(sr2a - c(1000, 1000, 0)) <= c(1, 2, 0)))
  
  ## may be off by small amount
  target.stats <- c(1500)
  nws22a <- san(nw ~ edges,
                target.stats = target.stats,
                constraints = ~bd(maxout = 3),
                control = control.san(SAN.maxit = 1, SAN.nsteps = 6e3))

  sr22a <- summary(nws22a ~ edges + degree(3) + degrange(4))
  expect_true(all(abs(sr22a - c(1500, 1000, 0)) <= c(2, 2, 0)))

  nw %v% "sex" <- rep(c("A", "B"), 500)
  target.stats <- c(500)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 2e3))

  sr3 <- summary(nws ~ edges + concurrent + nodematch("sex"))
  expect_equal(unname(sr3), c(500, 0, 0))

  target.stats <- c(500)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 2e3))

  sr4 <- summary(nws ~ edges + concurrent + nodematch("sex"))
  expect_equal(unname(sr4), c(500, 0, 500))
  
  ## may be off by small amount
  target.stats <- c(1500)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 6e3))

  sr5 <- summary(nws ~ edges + degree(3) + degrange(4) + nodematch("sex"))
  expect_true(all(abs(sr5 - c(1500, 1000, 0, 0)) <= c(2, 4, 0, 2)))
  
  ## may be off by small amount
  target.stats <- c(1500)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 6e3))

  sr6 <- summary(nws ~ edges + degree(3) + degrange(4) + nodematch("sex"))
  expect_true(all(abs(sr6 - c(1500, 1000, 0, 1500)) <= c(4, 4, 0, 4)))

  target.stats <- c(1000)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 4e3))

  sr7 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))
  expect_equal(unname(sr7), c(1000, 0, 0))

  target.stats <- c(1000)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 4e3))

  sr8 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))
  expect_equal(unname(sr8), c(1000, 0, 1000))
})

test_that("BDStratTNT works with bipartite networks and bd/blocks only", {
  nw <- network.initialize(900, bipartite = 100, directed = FALSE)

  target.stats <- c(100)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 4e2))

  sr <- summary(nws ~ edges + b1degree(1) + degrange(2))
  expect_equal(unname(sr), c(100, 100, 0))

  target.stats <- c(200)
  nws2 <- san(nws ~ edges,
              target.stats = target.stats,
              constraints = ~bd(maxout = 2),
              control = control.san(SAN.maxit = 1, SAN.nsteps = 4e2))

  sr2 <- summary(nws2 ~ edges + b1degree(2) + degrange(3))
  expect_true(all(abs(sr2 - c(200, 100, 0)) <= c(0, 0, 0)))

  target.stats <- c(300)
  nws22 <- san(nws2 ~ edges,
               target.stats = target.stats,
               constraints = ~bd(maxout = 3),
               control = control.san(SAN.maxit = 1, SAN.nsteps = 4e2))

  sr22 <- summary(nws22 ~ edges + b1degree(3) + degrange(4))
  expect_true(all(abs(sr22 - c(300, 100, 0)) <= c(0, 0, 0)))

  target.stats <- c(200)
  nws2a <- san(nw ~ edges,
               target.stats = target.stats,
               constraints = ~bd(maxout = 2),
               control = control.san(SAN.maxit = 1, SAN.nsteps = 8e2))

  sr2a <- summary(nws2a ~ edges + b1degree(2) + degrange(3))
  expect_true(all(abs(sr2a - c(200, 100, 0)) <= c(0, 0, 0)))

  target.stats <- c(300)
  nws22a <- san(nw ~ edges,
                target.stats = target.stats,
                constraints = ~bd(maxout = 3),
                control = control.san(SAN.maxit = 1, SAN.nsteps = 1.2e3))

  sr22a <- summary(nws22a ~ edges + b1degree(3) + degrange(4))
  expect_true(all(abs(sr22a - c(300, 100, 0)) <= c(0, 0, 0)))

  nw %v% "sex" <- c(rep(c("A", "B"), 50), rep(c("A", "B"), 450))
  target.stats <- c(100)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 4e2))

  sr3 <- summary(nws ~ edges + concurrent + nodematch("sex"))
  expect_equal(unname(sr3), c(100, 0, 0))

  target.stats <- c(100)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 4e2))

  sr4 <- summary(nws ~ edges + concurrent + nodematch("sex"))
  expect_equal(unname(sr4), c(100, 0, 100))

  target.stats <- c(300)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 1.2e3))

  sr5 <- summary(nws ~ edges + b1degree(3) + degrange(4) + nodematch("sex"))
  expect_true(all(abs(sr5 - c(300, 100, 0, 0)) <= c(0, 0, 0, 0)))

  target.stats <- c(300)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 1.2e3))

  sr6 <- summary(nws ~ edges + b1degree(3) + degrange(4) + nodematch("sex"))
  expect_true(all(abs(sr6 - c(300, 100, 0, 300)) <= c(0, 0, 0, 0)))

  target.stats <- c(200)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 8e2))

  sr7 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))
  expect_equal(unname(sr7), c(200, 0, 0))

  target.stats <- c(200)
  nws <- san(nw ~ edges,
             target.stats = target.stats,
             constraints = ~bd(maxout = 3) + blocks(attr = "sex", levels2 = !diag(TRUE, 2)),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 8e2))

  sr8 <- summary(nws ~ edges + degrange(4) + nodematch("sex"))
  expect_equal(unname(sr8), c(200, 0, 200))
})

test_that("BDStratTNT works with undirected unipartite networks", {
  nw <- network.initialize(1000, dir=FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W", 960))

  pmat <- matrix(1, 3, 3)
  diag(pmat) <- c(2, 2, 30)

  target.stats <- c(1000, 50, 50, 800)
  nws <- san(nw ~ edges + nodematch("race", levels = NULL, diff = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  
  # to test initialization code, redo the SAN run with different targets,
  # starting from the previous network
  pmat <- matrix(10, 3, 3)
  diag(pmat) <- c(7, 7, 20)

  target.stats <- c(1000, 125, 125, 350)
  nws2 <- san(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~strat(attr = "race", pmat = pmat))

  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  ## redo above with lower target stats and a nontrivial upper bound on degree
  nw <- network.initialize(1000, directed = FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W", 960))
  pmat <- matrix(1, 3, 3)
  diag(pmat) <- c(2, 2, 10)

  target.stats <- c(160, 20, 20, 100)
  nws <- san(nw ~ edges + nodematch("race", levels = NULL, diff = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5) + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)

  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10, 3, 3)
  diag(pmat) <- c(7, 7, 20)

  target.stats <- c(530, 30, 30, 450)
  nws2 <- san(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5) + strat(attr = "race", pmat = pmat))

  sr <- summary(nws2 ~ edges + nodematch("race", levels = NULL, diff = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)

  ## again but now also with a levels2 argument
  nw <- network.initialize(1000, directed = FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W", 960))
  nw %v% "sex" <- rep(c("M", "F"), 500)
  
  pmat <- matrix(1, 3, 3)
  diag(pmat) <- c(2, 2, 10)

  target.stats <- c(80, 10, 10, 50)
  nws <- san(nw ~ edges + nodematch("race", levels = NULL, diff = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5)
                            + blocks(attr = "sex", levels2 = diag(TRUE, 2))
                            + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)

  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10, 3, 3)
  diag(pmat) <- c(7, 7, 20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5)
                             + blocks(attr = "sex", levels2 = diag(TRUE, 2))
                             + strat(attr = "race", pmat = pmat))

  sr <- summary(nws2 ~ edges + nodematch("race", levels = NULL, diff = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), 0)

  ## this time with only same-sex ties
  nw <- network.initialize(1000, directed = FALSE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W", 960))
  nw %v% "sex" <- rep(c("M", "F"), 500)

  pmat <- matrix(1, 3, 3)
  diag(pmat) <- c(2, 2, 10)

  target.stats <- c(80, 10, 10, 50)
  nws <- san(nw ~ edges + nodematch("race", levels = NULL, diff = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5)
                            + blocks(attr = "sex", levels2 = !diag(TRUE, 2))
                            + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))

  # to test initialization code, redo the SAN run with different targets, starting from the previous network
  pmat <- matrix(10, 3, 3)
  diag(pmat) <- c(7, 7, 20)

  target.stats <- c(285, 15, 15, 230)
  nws2 <- san(nws ~ edges + nodematch("race", levels = NULL, diff = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5)
                             + blocks(attr = "sex", levels2 = !diag(TRUE, 2))
                             + strat(attr = "race", pmat = pmat))

  sr <- summary(nws2 ~ edges + nodematch("race",levels=NULL, diff=TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), network.edgecount(nws2))
})

test_that("BDStratTNT works with bipartite networks", {
  nw <- network.initialize(900, bipartite = 100, directed = FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W", 860))

  pmat <- matrix(c(0, 100, 2, 2, 0, 2, 100, 100, 0), 3, 3, byrow = TRUE)

  target.stats <- c(0, 2, 100, 100, 0, 100, 2, 2, 0)
  nws <- san(nw ~ nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(100, 10, 0, 0, 100, 100, 10, 10, 0), 3, 3, byrow = TRUE)
  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)  
  nws2 <- san(nws ~ nodemix("race", levels2 = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~strat(attr = "race", pmat = pmat3))

  sr <- summary(nws2 ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  # redo above but with a nontrivial upper bound on degree
  nw <- network.initialize(900, bipartite = 100, directed = FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W", 860))

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0), 3, 3, byrow = TRUE)

  target.stats <- c(0, 2, 90, 45, 0, 45, 2, 2, 0)
  nws <- san(nw ~ nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5) + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ nodemix("race",levels2=TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0), 3, 3, byrow = TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)
  nws2 <- san(nws ~ nodemix("race", levels2 = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5) + strat(attr = "race", pmat = pmat3))

  sr <- summary(nws2 ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)

  # ditto but also with levels2
  nw <- network.initialize(900, bipartite = 100, directed = FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W", 860))
  nw %v% "sex" <- rep(c("M", "F"), 500)

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0), 3, 3, byrow = TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5)
                            + blocks(attr = "sex", levels2 = diag(TRUE, 2))
                            + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), 0)

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0), 3, 3, byrow = TRUE)

  pmat3 <- (pmat + pmat2)/2
  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race", levels2 = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5)
                             + blocks(attr = "sex", levels2 = diag(TRUE, 2))
                             + strat(attr = "race", pmat = pmat3))

  sr <- summary(nws2 ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), 0)
  
  # this time with only same-sex ties
  nw <- network.initialize(900, bipartite = 100, directed = FALSE)

  nw %v% "race" <- c(rep("B", 20), rep("W", 60), rep("A", 40), rep("B", 20), rep("W", 860))
  nw %v% "sex" <- rep(c("M", "F"), 500)

  pmat <- matrix(c(0, 45, 2, 2, 0, 2, 90, 45, 0), 3, 3, byrow = TRUE)

  target.stats <- round(c(0, 2, 90, 45, 0, 45, 2, 2, 0)/2)
  nws <- san(nw ~ nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             control = control.san(SAN.maxit = 1, SAN.nsteps = 5e3),
             constraints = ~bd(maxout = 5)
                            + blocks(attr = "sex", levels2 = !diag(TRUE, 2))
                            + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws ~ degrange(6))), 0)
  expect_equal(unname(summary(nws ~ nodematch("sex"))), network.edgecount(nws))

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(85, 10, 0, 0, 42, 43, 10, 10, 0), 3, 3, byrow = TRUE)
  pmat3 <- (pmat + pmat2)/2
  target.stats <- round(c(pmat2)/2)
  nws2 <- san(nws ~ nodemix("race", levels2 = TRUE),
              target.stats = target.stats,
              control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4),
              constraints = ~bd(maxout = 5)
                             + blocks(attr = "sex", levels2 = !diag(TRUE, 2))
                             + strat(attr = "race", pmat = pmat3))

  sr <- summary(nws2 ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
  expect_equal(unname(summary(nws2 ~ degrange(6))), 0)  
  expect_equal(unname(summary(nws2 ~ nodematch("sex"))), network.edgecount(nws2))
})

test_that("BDStratTNT works with impossible targets", {
  nw <- network.initialize(1000, directed = FALSE)

  nw %v% "race" <- c(rep("A", 30), rep("B", 30), rep("W", 940))
  nw %v% "sex"  <- rep(c("W", "X", "Y", "Z"), 250)

  levels2 <- matrix(0, 4, 4)
  levels2[1, 3] <- levels2[3, 1] <- levels2[2, 2] <- levels2[3, 4] <- levels2[4, 3] <- levels2[4, 4] <- 1
  levels2 <- levels2 > 0

  pmat <- matrix(c(25, 50, 5, 50, 25, 5, 5, 5, 100), 3, 3, byrow = TRUE)

  # impossible to hit these exactly
  target.stats <- c(211, 25, 50, 25, 5, 5, 100)
  nws <- san(nw ~ edges + nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             constraints = ~bd(maxout = 4)
                            + blocks(attr = "sex", levels2 = levels2)
                            + strat(attr = "race", pmat = pmat))

  sr <- summary(nws ~ edges + nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats + 1))
  expect_equal(unname(summary(nws ~ degrange(5))), 0)
  # and check sex nodemix
  srs <- summary(nws ~ nodemix("sex", levels2 = TRUE))
  expect_true(all(srs[as.logical(levels2[upper.tri(levels2, diag = TRUE)])] == 0))
})

test_that("BDStratTNT simulates reasonably", {
  for(deg_bound in c(1, 3)) {
    net_size <- 2000L

    nw <- network.initialize(net_size, directed = FALSE)

    vattr <- sample(c("A", "B", "C"), net_size, TRUE)
    sex <- sample(c("X", "Y", "Z"), net_size, TRUE)

    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex

    levels2 <- matrix(c(1, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3)
    levels2 <- levels2 > 0

    pmat <- 1 - matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3)

    nw_sim <- nw

    for(i in 1:2) {
      nw_sim <- simulate(nw_sim ~ edges,
                         coef = c(0),
                         constraints = ~bd(maxout = deg_bound)
                                        + blocks(attr = "sex", levels2 = levels2)
                                        + strat(attr = "vattr", pmat = pmat),
                         output = "network")
      summ_stats <- summary(nw_sim ~ nodemix("vattr", levels2 = TRUE)
                                     + nodemix("sex", levels2 = TRUE)
                                     + degrange(deg_bound + 1))

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
  for(deg_bound in c(1, 3)) {
    net_size <- 2000L

    nw <- network.initialize(net_size, directed = FALSE)

    vattr <- sample(c("A", "B", "C"), net_size, TRUE)
    sex <- sample(c(1, 2, 3), net_size, TRUE)

    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE

    nw %v% "vattr" <- vattr
    nw %v% "sex" <- sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)

    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(1, 2, 2, 4), c(5, 2, 3, 4))] <- TRUE
    blocks_levels_2 <- blocks_levels_2 | t(blocks_levels_2)

    levels2 <- matrix(c(1, 0, 1, 0, 0, 0, 1, 0, 0), 3, 3)
    levels2 <- levels2 > 0

    pmat <- 1 - matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 0), 3, 3)

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

    for(i in 1:2) {
      nw_sim <- simulate(nw_sim ~ edges,
                         coef = c(0),
                         constraints = ~bd(attribs = attribs, maxout = maxout)
                                        + blocks(~blocks_attr, levels2 = blocks_levels_2)
                                        + strat(attr = "vattr", pmat = pmat),
                         output = "network")

      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr", levels2 = TRUE))
      expect_true(all(summ_stats_vattr[c(1,3)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,3)] > 0))

      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr", levels2 = TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3,5,10,11)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3,5,10,11)] > 0))

      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                    to = factor(sex[c(el[, c(2, 1)])], levels = seq_len(3)))
      expect_true(all(degs <= maxout))
    }
  }
})

test_that("BDStratTNT simulates reasonably with bipartite heterogeneous degree bounds", {
  for(deg_bound in c(1, 3)) {
    net_size <- 2000L
    bip <- 700L

    nw <- network.initialize(net_size, directed = FALSE, bipartite = bip)

    vattr <- c(sample(c("A", "B", "C", "D"), bip, TRUE), sample(c("X", "Y", "Z"), net_size - bip, TRUE))
    sex <- c(sample(c(1, 2, 3, 4, 5), bip, TRUE), sample(c(6, 7, 8, 9, 10, 11), net_size - bip, TRUE))

    attribs <- matrix(FALSE, nrow = net_size, ncol = length(unique(sex)))
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    

    nw %v% "vattr" <- vattr
    nw %v% "sex" <- sex
    nw %v% "blocks_attr" <- c(sample(c(1, 2, 3), bip, TRUE), sample(c(4, 5, 6, 7), net_size - bip, TRUE))

    blocks_levels_2 <- matrix(FALSE, nrow = 3, 4)
    blocks_levels_2[cbind(c(3, 2, 2), c(1, 2, 3))] <- TRUE

    levels2 <- matrix(as.logical(round(runif(11*11))), nrow = 11, ncol = 11)
    levels2 <- levels2 | t(levels2)

    pmat <- 1 - matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1), nrow = 4, ncol = 3)

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

    for(i in 1:2) {
      nw_sim <- simulate(nw_sim ~ edges,
                         coef = c(0),
                         constraints = ~bd(attribs = attribs, maxout = maxout)
                                        + blocks(~blocks_attr, levels2 = blocks_levels_2)
                                        + strat(attr = "vattr", pmat = pmat),
                         output = "network")

      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr", levels2 = TRUE))
      expect_true(all(summ_stats_vattr[c(1, 5, 7, 12)] == 0))
      expect_true(all(summ_stats_vattr[-c(1, 5, 7, 12)] > 0))

      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr", levels2 = TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3, 5, 8)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3, 5, 8)] > 0))

      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                    to = factor(sex[c(el[, c(2, 1)])], levels = seq_len(11)))
      expect_true(all(degs <= maxout))
    }
  }
})

test_that("BDStratTNT simulates reasonably with directed heterogeneous degree bounds", {
  for(deg_bound in c(1, 3)) {
    net_size <- 2000L

    nw <- network.initialize(net_size, directed = TRUE)

    vattr <- sample(c("A", "B", "C"), net_size, TRUE)
    sex <- sample(c(1, 2, 3), net_size, TRUE)

    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE

    nw %v% "vattr" <- vattr
    nw %v% "sex" <- sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)

    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(5, 2, 2, 4), c(1, 2, 3, 4))] <- TRUE

    levels2 <- matrix(c(1, 0, 0, 0, 0, 1, 1, 0, 0), 3, 3)
    levels2 <- levels2 > 0

    pmat <- 1 - matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 0), 3, 3)

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

    for(i in 1:2) {
      nw_sim <- simulate(nw_sim ~ edges,
                         coef = c(0),
                         constraints = ~bd(attribs = attribs, maxout = maxout, maxin = maxin)
                                        + blocks(~blocks_attr, levels2 = blocks_levels_2)
                                        + strat(attr = "vattr", pmat = pmat),
                         output = "network")

      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr", levels2 = TRUE))
      expect_true(all(summ_stats_vattr[c(1, 8)] == 0))
      expect_true(all(summ_stats_vattr[-c(1, 8)] > 0))

      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr" , levels2 = TRUE))
      expect_true(all(summ_stats_blocks_attr[c(5, 8, 14, 22)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(5, 8, 14, 22)] > 0))

      el <- as.edgelist(nw_sim)
      out_degs <- table(from = factor(c(el[, 1]), levels = seq_len(net_size)),
                        to = factor(sex[c(el[, 2])], levels = seq_len(3)))
      in_degs <- table(from = factor(c(el[, 2]), levels = seq_len(net_size)),
                       to = factor(sex[c(el[, 1])], levels = seq_len(3)))
      expect_true(all(out_degs <= maxout))
      expect_true(all(in_degs <= maxin))
    }
  }
})

test_that("BDStratTNT works with degree bound saturation", {
  nw <- network.initialize(900, directed = FALSE)

  nw %v% "race" <- rep(c("A", "B", "C"), times = c(30, 30, 840))
  nw %v% "sex" <- rep(c("X", "Y", "Z"), length.out = 900)

  pmat <- matrix(1, 3, 3)

  # mix.race.A.A mix.race.A.B mix.race.B.B mix.race.A.C mix.race.B.C mix.race.C.C
  target.stats <- c(425, 10, 5, 10, 0, 0, 400.01)

  nws <- san(nw ~ edges + nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             constraints = ~bd(maxout = 1)
                            + blocks(attr = "sex", levels2 = matrix(c(TRUE, FALSE, FALSE, FALSE, FALSE,
                                                                      TRUE, FALSE, TRUE, FALSE), 3, 3))
                            + strat(attr = "race", pmat = pmat),
             control = control.san(SAN.invcov.diag = TRUE, SAN.maxit = 4, SAN.nsteps = 5e4))
  sr <- summary(nws ~ edges + nodemix("race", levels2 = TRUE))

  expect_true(all(abs(sr - target.stats) <= pmax(1, 0.05*target.stats)))
  expect_true(all(summary(nws ~ concurrent + nodemix("sex", levels2=c(1, 5))) == 0))
})

test_that("BDStratTNT constrains undirected appropriately", {
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(c("A", "B", "C", "D", "E"), each = 20)
  nw %v% "strat_attr" <- rep(1:3, length.out = 100)
  nw[cbind(1:10, 30:21)] <- 1
  nw[cbind(44:53, 99:90)] <- 1
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~attr, levels2 = c(2, 13))
                                 + strat(~strat_attr, pmat = matrix(2 + runif(9), 3, 3)))

  expect_true(all(nws[cbind(1:10, 30:21)] == 1))
  expect_true(all(nws[cbind(44:53, 99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 13)) == 10)
  expect_true(summary(nws ~ edges) > 1000)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(~attr, levels2 = c(2, 13))
                                 + strat(~strat_attr, pmat = matrix(2 + runif(9), 3, 3)))

  expect_true(all(nws[cbind(1:10, 30:21)] == 1))
  expect_true(all(nws[cbind(44:53, 99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 13)) == 10)
  expect_true(summary(nws ~ edges) > 30)

  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(c("B", "A", "C", "D", "E"), each = 20)
  nw %v% "strat_attr" <- rep(1:3, length.out = 100)
  nw[cbind(1:10, 30:21)] <- 1
  nw[cbind(44:53, 99:90)] <- 1
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~attr, levels2 = c(2, 13))
                                 + strat(~strat_attr, pmat = matrix(2 + runif(9), 3, 3)))

  expect_true(all(nws[cbind(1:10, 30:21)] == 1))
  expect_true(all(nws[cbind(44:53, 99:90)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 2)) == 10)
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 13)) == 10)
  expect_true(summary(nws ~ edges) > 1000)
})

test_that("BDStratTNT constrains bipartite appropriately", {
  nw <- network.initialize(100, bipartite = 30, directed = FALSE)
  nw %v% "attr" <- c(rep(c("A", "B", "C"), each = 10),
                     rep(c("D", "E", "F", "G"), times = c(20, 20, 20, 10)))
  nw %v% "strat_attr" <- rep(1:3, length.out = 100)
  nw[cbind(1:10, 100:91)] <- 1
  nw[cbind(25:21, 60:56)] <- 1
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~attr, levels2 = c(6,10))
                                 + strat(~strat_attr, pmat = matrix(2 + runif(9), 3, 3)))

  expect_true(all(nws[cbind(1:10, 100:91)] == 1))
  expect_true(all(nws[cbind(25:21, 60:56)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 6)) == 5)
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 10)) == 10)
  expect_true(summary(nws ~ edges) > 500)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(~attr, levels2 = c(6, 10))
                                 + strat(~strat_attr, pmat = matrix(2 + runif(9), 3, 3)))

  expect_true(all(nws[cbind(1:10, 100:91)] == 1))
  expect_true(all(nws[cbind(25:21, 60:56)] == 1))
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 6)) == 5)
  expect_true(summary(nws ~ nodemix(~attr, levels2 = 10)) == 10)
  expect_true(summary(nws ~ edges) > 20)
})

test_that("BDStratTNT handles undirected arguments correctly", {
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "bd_attr" <- rep(1:3, length.out = 100)
  nw %v% "strat_attr" <- rep(1:7, length.out = 100)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 8)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -1)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 8)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(attr = ~bd_attr,
                                        levels2 = matrix(c(FALSE, TRUE, FALSE, TRUE, rep(FALSE, 5)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -2)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(FALSE, TRUE, FALSE,
                                                                              TRUE, rep(FALSE, 5)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 8)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -1)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 8)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(FALSE, TRUE, FALSE,
                                                                              TRUE, rep(FALSE, 5)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -2)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(FALSE, TRUE, FALSE, TRUE, rep(FALSE, 5)), 3, 3)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
})

test_that("BDStratTNT handles bipartite arguments correctly", {
  nw <- network.initialize(100, directed = FALSE, bipartite = 30)
  nw %v% "bd_attr" <- c(rep(1:3, length.out = 30), rep(6:10, length.out = 70))
  nw %v% "strat_attr" <- rep(1:7, length.out = 100)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(attr = ~bd_attr,
                                        levels2 = matrix(c(TRUE, rep(FALSE, 14)), nrow = 3,ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -1)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(attr = ~bd_attr,
                                          levels2 = matrix(c(TRUE, rep(FALSE, 14)), nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(attr = ~bd_attr,
                                        levels2 = matrix(c(FALSE, TRUE, FALSE, FALSE,
                                                           rep(FALSE, 11)), nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2=2)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2=-2)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = 1)
                                 + blocks(attr = ~bd_attr,
                                          levels2 = matrix(c(FALSE, TRUE, FALSE, FALSE,
                                                             rep(FALSE, 11)), nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 14)),
                                                                            nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -1)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + bd(maxout = 1)
                                 + blocks(attr = ~bd_attr, levels2 = matrix(c(TRUE, rep(FALSE, 14)),
                                                                            nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 1)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + blocks(attr = ~bd_attr,
                                          levels2 = matrix(c(FALSE, TRUE, FALSE, FALSE,
                                                             rep(FALSE, 11)), nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -2)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~strat(attr = "strat_attr", pmat = matrix(2 + runif(7*7), 7, 7))
                                 + bd(maxout = 1)
                                 + blocks(attr = ~bd_attr,
                                          levels2 = matrix(c(FALSE, TRUE, FALSE, FALSE,
                                                             rep(FALSE, 11)), nrow = 3, ncol = 5)),
                  control = list(MCMC.prop.weights = "BDStratTNT"))
  expect_true(summary(nws ~ nodemix(~bd_attr, levels2 = 2)) == 0)
  expect_true(all(summary(nws ~ nodefactor(~bd_attr, levels = TRUE)) > 0))
  expect_true(summary(nws ~ concurrent) == 0)
})

test_that("BDStratTNT handles atypical levels specifications correctly", {
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "bd_attr" <- rep(1:3, length.out = 100)
  nw %v% "strat_attr" <- rep(1:5, length.out = 100)
  pmat <- matrix(2 + runif(25), 5, 5)

  ## should be unconstrained
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, levels = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  ## should also be unconstrained
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, levels = FALSE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  ## any pairing with a 3 should be allowed, with all other pairings forbidden
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, levels = I(c(1, 2, 4, 6)), levels2 = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = c(4, 5, 6))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -c(4, 5, 6))) == 0))

  ## only 2-2 pairings should be allowed
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, levels = I(c(1, 2, 3, 4, 6)), levels2 = -3)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = c(3))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -c(3))) == 0))

  ## similar bipartite tests
  nw <- network.initialize(100, directed = FALSE, bipartite = 30)
  nw %v% "bd_attr" <- c(rep(1:3, length.out = 30), rep(10:16, length.out = 70))
  nw %v% "strat_attr" <- c(rep(1:5, length.out = 30), rep(1:4, length.out = 70))
  pmat <- matrix(2 + runif(20), nrow = 5, ncol = 4)

  ## should be unconstrained
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = TRUE, b2levels = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  ## should also be unconstrained
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = FALSE, b2levels = FALSE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = FALSE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b2levels = FALSE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = FALSE, b2levels = FALSE, levels2 = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = FALSE, levels2 = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b2levels = FALSE, levels2 = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = TRUE)) > 0))

  ## any pairing with a 3 should be allowed, with all other pairings forbidden
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = I(c(1, 2, 4, 6)), levels2 = TRUE)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = 3*(1:7))) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -3*(1:7))) == 0))

  ## only 1-14 pairings should be allowed
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~blocks(~bd_attr, b1levels = I(c(1, 2, 3, 4, 6)), levels2 = -21)
                                 + strat(attr = ~strat_attr, pmat = pmat))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, b1levels = I(1), b2levels = I(14), levels2 = TRUE)) > 0))
  expect_true(all(summary(nws ~ nodemix(~bd_attr, levels2 = -c(13))) == 0))
})

test_that("BDStratTNT works with directed networks", {
  nw <- network.initialize(1000, directed = TRUE)

  nw %v% "race" <- c(rep("A", 20), rep("B", 20), rep("W", 960))

  pmat <- matrix(c(100, 350, 0, 10, 100, 0, 100, 0, 840), 3, 3, byrow = TRUE)

  target.stats <- c(100, 10, 100, 350, 100, 0, 0, 0, 840)
  nws <- san(nw ~ nodemix("race", levels2 = TRUE),
             target.stats = target.stats,
             constraints = ~bd(maxout = 40, maxin = 40)
                            + strat(pmat = pmat, attr = "race"),
             control = control.san(SAN.maxit = 1, SAN.nsteps = 1e4))

  sr <- summary(nws ~ nodemix("race", levels2 = TRUE))
  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))

  # redo with different targets, starting from previous network
  pmat2 <- matrix(c(50, 50, 350, 50, 50, 100, 50, 400, 400), 3, 3, byrow = TRUE)

  pmat3 <- (pmat + pmat2)/2

  target.stats <- c(pmat2)
  nws2 <- san(nws ~ nodemix("race", levels2 = TRUE),
              target.stats = target.stats,
              constraints = ~bd(maxout = 40, maxin = 40)
                             + strat(pmat = pmat3, attr = "race"),
              control = control.san(SAN.maxit = 1, SAN.nsteps = 2e4))
  sr <- summary(nws2 ~ nodemix("race", levels2 = TRUE))

  expect_true(all(abs(sr - target.stats) <= 0.05*target.stats))
})

test_that("BDStratTNT simulates directed reasonably", {
  net_size <- 1000L

  nw <- network.initialize(net_size, directed = TRUE)

  vattr <- sample(c("A", "B", "C"), net_size, TRUE)

  nw %v% "vattr" <- vattr
  nw %v% "sex" <- sample(c("X", "Y", "Z"), net_size, TRUE)

  pmat <- 1 - matrix(c(1, 0, 0, 1, 1, 0, 0, 1, 0), 3, 3)

  nw_sim <- nw

  for(i in c(1, 3)) {
    nw_sim <- simulate(nw_sim ~ edges,
                       coef = c(0),
                       constraints = ~bd(maxout = i, maxin = i + 1)
                                      + blocks(attr = "sex",
                                               levels2 = matrix(c(TRUE, FALSE, TRUE, FALSE, FALSE,
                                                                  TRUE, FALSE, TRUE, FALSE), 3, 3))
                                      + strat(attr = "vattr", pmat = pmat),
                       output = "network")
    summ_stats <- summary(nw_sim ~ nodemix("vattr", levels2 = TRUE)
                                   + nodemix("sex", levels2 = TRUE)
                                   + odegrange(i + 1)
                                   + idegrange(i + 2))

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

test_that("BDStratTNT handles undirected heterogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = FALSE)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(10*10), nrow = 10, ncol = 10)
  pmat <- pmat + t(pmat)

  levels2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- matrix(round(deg_bound*runif(net_size*7)), nrow = net_size)
  bd_attr <- matrix(FALSE, nrow = net_size, ncol = 7)
  bd_attr[cbind(seq_len(net_size), 1 + (seq_len(net_size) %% 7))] <- TRUE  
  bd_attr_flat <- rep(c(2:7,1), length.out = net_size)
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(attr = bd_attr, maxout = maxout)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                to = factor(bd_attr_flat[c(el[, c(2, 1)])], levels = seq_len(7)))
  expect_true(all(degs <= maxout))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(attr = bd_attr, maxout = maxout)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                to = factor(bd_attr_flat[c(el[, c(2, 1)])], levels = seq_len(7)))
  expect_true(all(degs <= maxout))
})

test_that("BDStratTNT handles directed heterogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = TRUE)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(10*10), nrow = 10, ncol = 10)

  levels2 <- matrix(c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- matrix(round(deg_bound*runif(net_size*7)), nrow = net_size)
  maxin <- matrix(round(deg_bound*runif(net_size*7)), nrow = net_size)
  bd_attr <- matrix(FALSE, nrow = net_size, ncol = 7)
  bd_attr[cbind(seq_len(net_size), 1 + (seq_len(net_size) %% 7))] <- TRUE  
  bd_attr_flat <- rep(c(2:7,1), length.out = net_size)
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(attr = bd_attr, maxout = maxout, maxin = maxin)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  out_degs <- table(from = factor(c(el[, 1]), levels = seq_len(net_size)),
                    to = factor(bd_attr_flat[c(el[, 2])], levels = seq_len(7)))
  expect_true(all(out_degs <= maxout))
  in_degs <- table(from = factor(c(el[, 2]), levels = seq_len(net_size)),
                   to = factor(bd_attr_flat[c(el[, 1])], levels = seq_len(7)))
  expect_true(all(in_degs <= maxin))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(attr = bd_attr, maxout = maxout, maxin = maxin)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  out_degs <- table(from = factor(c(el[, 1]), levels = seq_len(net_size)),
                    to = factor(bd_attr_flat[c(el[, 2])], levels = seq_len(7)))
  expect_true(all(out_degs <= maxout))
  in_degs <- table(from = factor(c(el[, 2]), levels = seq_len(net_size)),
                   to = factor(bd_attr_flat[c(el[, 1])], levels = seq_len(7)))
  expect_true(all(in_degs <= maxin))
})

test_that("BDStratTNT handles bipartite heterogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  bip_size <- 5
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = FALSE, bipartite = bip_size)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(5*10), nrow = 5, ncol = 10)

  levels2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- matrix(round(deg_bound*runif(net_size*7)), nrow = net_size)
  bd_attr <- matrix(FALSE, nrow = net_size, ncol = 7)
  bd_attr[cbind(seq_len(net_size), 1 + (seq_len(net_size) %% 7))] <- TRUE  
  bd_attr_flat <- rep(c(2:7,1), length.out = net_size)
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(attr = bd_attr, maxout = maxout)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                to = factor(bd_attr_flat[c(el[, c(2, 1)])], levels = seq_len(7)))
  expect_true(all(degs <= maxout))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(attr = bd_attr, maxout = maxout)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  degs <- table(from = factor(c(el), levels = seq_len(net_size)),
                to = factor(bd_attr_flat[c(el[, c(2, 1)])], levels = seq_len(7)))
  expect_true(all(degs <= maxout))
})

test_that("BDStratTNT handles undirected homogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = FALSE)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(10*10), nrow = 10, ncol = 10)
  pmat <- pmat + t(pmat)

  levels2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- deg_bound
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = maxout)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  degs <- tabulate(c(el), nbins = net_size)
  expect_true(all(degs <= maxout))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(maxout = maxout)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  degs <- tabulate(c(el), nbins = net_size)
  expect_true(all(degs <= maxout))
})

test_that("BDStratTNT handles directed homogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = TRUE)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(10*10), nrow = 10, ncol = 10)

  levels2 <- matrix(c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- deg_bound
  maxin <- deg_bound
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = maxout, maxin = maxin)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  out_degs <- tabulate(c(el[, 1]), nbins = net_size)
  expect_true(all(out_degs <= maxout))
  in_degs <- tabulate(c(el[, 2]), nbins = net_size)
  expect_true(all(in_degs <= maxin))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(maxout = maxout, maxin = maxin)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  out_degs <- tabulate(c(el[, 1]), nbins = net_size)
  expect_true(all(out_degs <= maxout))
  in_degs <- tabulate(c(el[, 2]), nbins = net_size)
  expect_true(all(in_degs <= maxin))
})

test_that("BDStratTNT handles bipartite homogeneous degree bound saturation correctly in simulation context", {
  net_size <- 20
  bip_size <- 5
  deg_bound <- 2
  nw <- network.initialize(net_size, directed = FALSE, bipartite = bip_size)
  nw %v% "strat_attr" <- rep(letters[1:10], length.out = net_size)
  nw %v% "blocks_attr" <- rep(1:3, length.out = net_size)

  pmat <- matrix(runif(5*10), nrow = 5, ncol = 10)

  levels2 <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE),
                    nrow = 3, byrow = TRUE)

  maxout <- deg_bound
  nws <- simulate(nw ~ edges,
                  coef = c(0),
                  constraints = ~bd(maxout = maxout)
                                 + blocks(attr = ~blocks_attr, levels2 = levels2)
                                 + strat(attr = ~strat_attr, pmat = pmat),
                  control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws)
  degs <- tabulate(c(el), nbins = net_size)
  expect_true(all(degs <= maxout))

  ## restart to test initialization
  nws2 <- simulate(nws ~ edges,
                   coef = c(0),
                   constraints = ~bd(maxout = maxout)
                                  + blocks(attr = ~blocks_attr, levels2 = levels2)
                                  + strat(attr = ~strat_attr, pmat = pmat),
                   control = list(MCMC.burnin = 1e4))
  ## check constraints
  expect_true(all(summary(nws2 ~ nodemix(~blocks_attr, levels2 = levels2)) == 0))
  el <- as.edgelist(nws2)
  degs <- tabulate(c(el), nbins = net_size)
  expect_true(all(degs <= maxout))
})

test_that("all free dyads vary with a blocks constraint", {
  net_size <- 17
  bip_size <- 9
  nsim <- 200

  blocks_levels_2 <- matrix(FALSE, nrow = 3, ncol = 3)
  blocks_levels_2[1,1] <- TRUE
  blocks_levels_2[2,3] <- TRUE

  net_attrs_list <- list(list(n = net_size, directed = FALSE, bipartite = FALSE),
                         list(n = net_size, directed = TRUE, bipartite = FALSE),
                         list(n = net_size, directed = FALSE, bipartite = bip_size))

  for (i in seq_along(net_attrs_list)) {
    ## construct network
    nw <- do.call(network.initialize, net_attrs_list[[i]])

    ## construct logical matrix of dyads in the network
    dyad_mat <- !diag(net_size)
    if (is.bipartite(nw)) {
      dyad_mat[-seq_len(bip_size),] <- FALSE
      dyad_mat[,seq_len(bip_size)] <- FALSE
    } else if (!is.directed(nw)) {
      dyad_mat[lower.tri(dyad_mat)] <- FALSE
    }

    ## check dyad count
    expect_equal(sum(dyad_mat), network.dyadcount(nw))

    ## construct logical vector of free dyads in the network
    free_dyads <- !diag(net_size)
    free_dyads[seq_len(net_size) %% 3 == 1, seq_len(net_size) %% 3 == 1] <- FALSE
    free_dyads[seq_len(net_size) %% 3 == 2, seq_len(net_size) %% 3 == 0] <- FALSE
    if (!is.directed(nw) && !is.bipartite(nw)) {
      free_dyads <- free_dyads & t(free_dyads)
    }
    free_dyads <- free_dyads[which(dyad_mat)]

    ## restrict size of dyad_mat in bipartite case
    if (is.bipartite(nw)) {
      dyad_mat <- dyad_mat[seq_len(bip_size), -seq_len(bip_size)]
    }

    ## perform simulation
    stats <- simulate(nw ~ edges,
                      coef = c(0),
                      monitor = ~nodemix(~seq_len(net_size), levels2 = dyad_mat),
                      nsim = nsim,
                      constraints = "BDStratTNT" ~ blocks(~rep(1:3, length.out = net_size),
                                                          levels2 = blocks_levels_2),
                      output = "stats",
                      control = list(MCMC.burnin = 1e3, MCMC.interval = 1e3))

    ## drop edges statistic
    stats <- stats[,-1]

    ## check number of stats
    expect_equal(NCOL(stats), network.dyadcount(nw))

    ## confirm that free dyads take two states and fixed dyads take one state
    for (i in seq_len(NCOL(stats))) {
      expect_true(length(unique(stats[, i])) == 1 + free_dyads[i])
    }
  }
})
