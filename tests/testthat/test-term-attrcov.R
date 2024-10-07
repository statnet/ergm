#  File tests/testthat/test-term-attrcov.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

test_that("attrcov works with undirected unipartite networks", {
  md <- matrix(runif(50*50), 50, 50)
  md[lower.tri(md)] <- t(md)[lower.tri(md)]
  mi <- matrix(sample(-10:10, 50*50, TRUE), 50, 50)
  mi[lower.tri(mi)] <- t(mi)[lower.tri(mi)]

  attr <- sample(rep(1:50, each=10))

  nw <- network.initialize(500, dir=FALSE)
  nw %v% "attr" <- attr
  
  nw_san <- san(nw ~ edges, target.stats=500)

  el <- as.edgelist(nw_san)
  el_attr <- matrix(attr[el], ncol=2)
  
  expect_true(sum(md[el_attr]) == summary(nw_san ~ attrcov(~attr, md)))
  expect_true(sum(mi[el_attr]) == summary(nw_san ~ attrcov(~attr, mi)))
  
  wto <- sample(seq_len(NROW(el)), floor(NROW(el)/2))
  
  changes <- rbind(el, el[wto,])
  changes <- changes[sample(seq_len(NROW(changes))),]
  
  rvd <- ergm.godfather(nw ~ attrcov(~attr, md), changes=changes)
  expect_true(sum(md[el_attr[-wto,]]) == rvd)

  rvi <- ergm.godfather(nw ~ attrcov(~attr, mi), changes=changes)
  expect_true(sum(mi[el_attr[-wto,]]) == rvi)  

  mw <- matrix(runif(40*50), 40, 50)
  expect_error(summary(nw ~ attrcov(~attr, mw)))

  mw2 <- matrix(runif(50*20), 50, 20)
  expect_error(summary(nw ~ attrcov(~attr, mw2)))

  mw3 <- matrix(runif(500*5), 500, 5)
  expect_error(summary(nw ~ attrcov(~attr, mw3)))

  mw4 <- matrix(runif(100*100), 100, 100)
  mw4[lower.tri(mw4)] <- t(mw4)[lower.tri(mw4)]
  expect_error(summary(nw ~ attrcov(~attr, mw4)))
})

test_that("attrcov works with undirected bipartite networks", {
  md <- matrix(runif(30*50), 30, 50)
  mi <- matrix(sample(-10:10, 30*50, TRUE), 30, 50)

  attr <- c(sample(rep(1:30, length.out=200)), sample(rep(1:50, length.out=300)))

  nw <- network.initialize(500, dir=FALSE, bip=200)
  nw %v% "attr" <- attr
  
  nw_san <- san(nw ~ edges, target.stats=500)

  el <- as.edgelist(nw_san)
  el_attr <- matrix(attr[el], ncol=2)
  
  expect_true(sum(md[el_attr]) == summary(nw_san ~ attrcov(~attr, md)))
  expect_true(sum(mi[el_attr]) == summary(nw_san ~ attrcov(~attr, mi)))
  
  wto <- sample(seq_len(NROW(el)), floor(NROW(el)/2))
  
  changes <- rbind(el, el[wto,])
  changes <- changes[sample(seq_len(NROW(changes))),]
  
  rvd <- ergm.godfather(nw ~ attrcov(~attr, md), changes=changes)
  expect_true(sum(md[el_attr[-wto,]]) == rvd)

  rvi <- ergm.godfather(nw ~ attrcov(~attr, mi), changes=changes)
  expect_true(sum(mi[el_attr[-wto,]]) == rvi)  
  
  mw <- matrix(runif(40*50), 40, 50)
  expect_error(summary(nw ~ attrcov(~attr, mw)))  

  mw2 <- matrix(runif(30*40), 30, 40)
  expect_error(summary(nw ~ attrcov(~attr, mw2)))

  mw3 <- matrix(runif(300*5), 300, 5)
  expect_error(summary(nw ~ attrcov(~attr, mw3)))
  
  mw4 <- matrix(runif(100*100), 100, 100)
  expect_error(summary(nw ~ attrcov(~attr, mw4)))  
})

test_that("attrcov works with directed networks", {
  md <- matrix(runif(50*50), 50, 50)
  mi <- matrix(sample(-10:10, 50*50, TRUE), 50, 50)

  attr <- sample(rep(1:50, each=10))

  nw <- network.initialize(500, dir=TRUE)
  nw %v% "attr" <- attr
  
  nw_san <- san(nw ~ edges, target.stats=500)

  el <- as.edgelist(nw_san)
  el_attr <- matrix(attr[el], ncol=2)
  
  expect_true(sum(md[el_attr]) == summary(nw_san ~ attrcov(~attr, md)))
  expect_true(sum(mi[el_attr]) == summary(nw_san ~ attrcov(~attr, mi)))
  
  wto <- sample(seq_len(NROW(el)), floor(NROW(el)/2))
  
  changes <- rbind(el, el[wto,])
  changes <- changes[sample(seq_len(NROW(changes))),]
  
  rvd <- ergm.godfather(nw ~ attrcov(~attr, md), changes=changes)
  expect_true(sum(md[el_attr[-wto,]]) == rvd)

  rvi <- ergm.godfather(nw ~ attrcov(~attr, mi), changes=changes)
  expect_true(sum(mi[el_attr[-wto,]]) == rvi)  

  mw <- matrix(runif(40*50), 40, 50)
  expect_error(summary(nw ~ attrcov(~attr, mw)))

  mw2 <- matrix(runif(50*20), 50, 20)
  expect_error(summary(nw ~ attrcov(~attr, mw2)))

  mw3 <- matrix(runif(500*5), 500, 5)
  expect_error(summary(nw ~ attrcov(~attr, mw3)))

  mw4 <- matrix(runif(100*100), 100, 100)
  expect_error(summary(nw ~ attrcov(~attr, mw4)))
})
