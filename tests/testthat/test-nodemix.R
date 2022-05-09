#  File tests/testthat/test-nodemix.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
library(ergm)


# Test undirected network
data(faux.mesa.high)
m <- matrix(c(75, 0, 0, 1, 1, 1, 0, 33, 2, 4, 2, 1,
 0, 2, 23, 7, 6, 4, 1, 4, 7, 9, 1, 5, 1, 2, 6, 1, 17, 5,
 1, 1, 4, 5, 5, 6), 6, 6)  # Correct answer!
faux.mesa.high %e% "eattr" <- rep(1, network.edgecount(faux.mesa.high))

test_that("binary nodemix on undirected network faux.mesa.high", {
  expect_equal(m[upper.tri(m, diag=T)], summary(faux.mesa.high ~ nodemix("Grade", levels2=TRUE)), ignore_attr = TRUE)
})
test_that("weighted nodemix on undirected network faux.mesa.high", {
  expect_equal(m[upper.tri(m, diag=T)], summary(faux.mesa.high ~ nodemix("Grade", form="nonzero", levels2=TRUE), response="eattr"), ignore_attr = TRUE)
})

# directed network
data(sampson)
grpord<-c("Turks","Loyal","Outcasts")
m2 <- matrix(c(30, 9, 7, 5, 23, 1, 1, 2, 10), 3, 3,dimnames=list(From=grpord,To=grpord)) # Correct answer!
mixnames<-t(sapply(grpord,function(from) sapply(grpord,function(to) paste("mix.group",from,to,sep="."))))
samplike %e% "eattr" <- rep(1, network.edgecount(samplike))

test_that("binary nodemix failed test on directed network samplike", {
  expect_equal(c(m2), summary(samplike ~ nodemix("group", levels2=TRUE))[c(mixnames)], ignore_attr = TRUE)
})
test_that("weighted nodemix failed test on directed network samplike", {
  expect_equal(c(m2), summary(samplike ~ nodemix("group", form="nonzero", levels2=TRUE), response="eattr")[c(mixnames)], ignore_attr = TRUE)
})

# bipartite network
el <- cbind( c(17, 20, 22, 26, 19, 24, 16, 22, 18, 23, 28, 20,
               22, 23, 17, 21, 25, 21, 27, 16, 19, 18, 23),
           c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 10, 
             10, 11, 11))
mynw <- network(el, bipartite=15, directed=FALSE) 
mynw %v% "names" <- rep(letters[1:3], c(10,10,8))
m3 <- matrix(c(9, 1, 12, 1), 2, 2) # Correct answer!
mynw %e% "eattr" <- rep(1, network.edgecount(mynw))

test_that("binary nodemix failed test on bipartite network", {
  expect_equal(as.vector(m3), summary(mynw ~ nodemix("names", levels2=TRUE)), ignore_attr = TRUE)
})
test_that("weighted nodemix failed test on bipartite network", {
  expect_equal(as.vector(m3), summary(mynw ~ nodemix("names", form="nonzero", levels2=TRUE), response="eattr"), ignore_attr = TRUE)
})

net_size <- 1000
bip_size <- 300
attr_levels <- 6
b1attr_levels <- 5
b2attr_levels <- 7

## summary tests ##

test_that("nodemix undirected summary test", {

nw <- network.initialize(net_size, dir=FALSE)
nw <- san(nw ~ edges, target.stats=net_size)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- rep(seq_len(attr_levels), length.out=net_size)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(0L, attr_levels, attr_levels)
m_ind[upper.tri(m_ind, diag=TRUE)] <- seq_len(attr_levels*(attr_levels + 1)/2)
m_ind <- m_ind + t(m_ind) - diag(diag(m_ind))
levs <- -c(5,3)
levs2 <- c(1,15,6,3,7,4,20,18)
levs2withlevs <- c(8,2,7,6,4)

stats <- rep(0, attr_levels*(attr_levels + 1)/2)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
expect_equal(stats, summary(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr"), ignore_attr = TRUE)

stats <- rep(0, attr_levels*(attr_levels + 1)/2)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, summary(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr"), ignore_attr = TRUE)

expect_equal(stats, summary(nw ~ nodemix("vattr", levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs)), ignore_attr = TRUE)

})

## directed ##

test_that("nodemix directed summary test", {

nw <- network.initialize(net_size, dir=TRUE)
nw <- san(nw ~ edges, target.stats=net_size)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- rep(seq_len(attr_levels), length.out=net_size)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(seq_len(attr_levels*attr_levels), attr_levels, attr_levels)
levs <- -c(5,3)
levs2 <- c(1,15,33,3,27,24,30,18,5,9,11)
levs2withlevs <- c(8,12,7,6,14,4,3)

stats <- rep(0, attr_levels*attr_levels)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
expect_equal(stats, summary(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr"), ignore_attr = TRUE)

stats <- rep(0, attr_levels*attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, summary(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr"), ignore_attr = TRUE)

expect_equal(stats, summary(nw ~ nodemix("vattr", levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], summary(nw ~ nodemix("vattr", levels = levs, levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs)), ignore_attr = TRUE)

})

## bipartite ##

test_that("nodemix bipartite summary test", {

nw <- network.initialize(net_size, dir=FALSE, bip=bip_size)
nw <- san(nw ~ edges, target.stats=net_size)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- c(rep(seq_len(b1attr_levels), length.out=bip_size), rep(seq_len(b2attr_levels), length.out=net_size - bip_size))
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(seq_len(b1attr_levels*b2attr_levels), b1attr_levels, b2attr_levels)
b1levs <- c(1,2,4,5)
b2levs <- -c(3,6)
levs2 <- c(34,23,11,6,3,9,19,25,28,30)
levs2withblevs <- c(16,4,12,3,19,17,8)

stats <- rep(0, b1attr_levels*b2attr_levels)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
expect_equal(stats, summary(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="sum", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="sum"),response="eattr"), ignore_attr = TRUE)

stats <- rep(0, b1attr_levels*b2attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, summary(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="nonzero", levels2=TRUE),response="eattr"), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="nonzero"),response="eattr"), ignore_attr = TRUE)

expect_equal(stats, summary(nw ~ nodemix("vattr", levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[levs2], summary(nw ~ nodemix("vattr", levels2 = levs2)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2=TRUE)), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs)), ignore_attr = TRUE)

})

## godfather tests ##

## undirected ##

test_that("nodemix undirected godfather test", {

nw <- network.initialize(net_size, dir=FALSE)

nw <- san(nw ~ edges, target.stats=as.integer(3*net_size/2))
el <- as.edgelist(nw)

indices <- seq_len(NROW(el))

init_indices <- sample(indices, net_size, FALSE)

toggle_off <- sample(init_indices, as.integer(net_size/2), FALSE)

godfather_indices <- c(toggle_off, indices[-init_indices])

final_indices <- c(setdiff(init_indices, toggle_off), indices[-init_indices])

el_init <- el[sort(init_indices),,drop=FALSE]
attr(el_init, "n") <- net_size

nw <- network(el_init, type="edgelist", directed = FALSE)
wts <- matrix(runif(net_size*net_size), net_size, net_size)
wts <- (wts + t(wts))/2

eattr <- wts[el_init]
vattr <- rep(seq_len(attr_levels), length.out=net_size)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(0L, attr_levels, attr_levels)
m_ind[upper.tri(m_ind, diag=TRUE)] <- seq_len(attr_levels*(attr_levels + 1)/2)
m_ind <- m_ind + t(m_ind) - diag(diag(m_ind))
levs <- -c(5,3)
levs2 <- c(1,15,6,3,7,4,20,18)
levs2withlevs <- c(8,2,7,6,4)

el_off <- el[toggle_off,,drop=FALSE]
el_on <- el[indices[-init_indices],,drop=FALSE]
changes_godfather <- rbind(cbind(el_off, 0), cbind(el_on, wts[el_on]))
changes_godfather <- changes_godfather[sample(seq_len(NROW(changes_godfather))),,drop=FALSE]

el_final <- el[final_indices,,drop=FALSE]
el_final <- cbind(el_final, wts[el_final])
el_final <- el_final[sample(seq_len(NROW(el_final))),,drop=FALSE]

stats <- rep(0, attr_levels*(attr_levels + 1)/2)
types <- m_ind[matrix(vattr[el_final], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + wts[el_final[i,1:2,drop=FALSE]]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

stats <- rep(0, attr_levels*(attr_levels + 1)/2)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

changes_godfather <- changes_godfather[,1:2,drop=FALSE]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs),changes=changes_godfather), ignore_attr = TRUE)

})

## directed ##

test_that("nodemix directed godfather test", {

nw <- network.initialize(net_size, dir=TRUE)

nw <- san(nw ~ edges, target.stats=as.integer(3*net_size/2))
el <- as.edgelist(nw)

indices <- seq_len(NROW(el))

init_indices <- sample(indices, net_size, FALSE)

toggle_off <- sample(init_indices, as.integer(net_size/2), FALSE)

godfather_indices <- c(toggle_off, indices[-init_indices])

final_indices <- c(setdiff(init_indices, toggle_off), indices[-init_indices])

el_init <- el[sort(init_indices),,drop=FALSE]
attr(el_init, "n") <- net_size

nw <- network(el_init, type="edgelist", directed = TRUE)
wts <- matrix(runif(net_size*net_size), net_size, net_size)

eattr <- wts[el_init]
vattr <- rep(seq_len(attr_levels), length.out=net_size)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(seq_len(attr_levels*attr_levels), attr_levels, attr_levels)

levs <- -c(5,3)
levs2 <- c(1,15,33,3,27,24,30,18,5,9,11)
levs2withlevs <- c(8,12,7,6,14,4,3)

el_off <- el[toggle_off,,drop=FALSE]
el_on <- el[indices[-init_indices],,drop=FALSE]
changes_godfather <- rbind(cbind(el_off, 0), cbind(el_on, wts[el_on]))
changes_godfather <- changes_godfather[sample(seq_len(NROW(changes_godfather))),,drop=FALSE]

el_final <- el[final_indices,,drop=FALSE]
el_final <- cbind(el_final, wts[el_final])
el_final <- el_final[sample(seq_len(NROW(el_final))),,drop=FALSE]

stats <- rep(0, attr_levels*attr_levels)
types <- m_ind[matrix(vattr[el_final], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + wts[el_final[i,1:2,drop=FALSE]]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

stats <- rep(0, attr_levels*attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

changes_godfather <- changes_godfather[,1:2,drop=FALSE]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]], ergm.godfather(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs),changes=changes_godfather), ignore_attr = TRUE)

})

## bipartite ##

test_that("nodemix bipartite godfather test", {

nw <- network.initialize(net_size, dir=FALSE, bip=bip_size)

nw <- san(nw ~ edges, target.stats=as.integer(3*net_size/2))
el <- as.edgelist(nw)

indices <- seq_len(NROW(el))

init_indices <- sample(indices, net_size, FALSE)

toggle_off <- sample(init_indices, as.integer(net_size/2), FALSE)

godfather_indices <- c(toggle_off, indices[-init_indices])

final_indices <- c(setdiff(init_indices, toggle_off), indices[-init_indices])

el_init <- el[sort(init_indices),,drop=FALSE]
attr(el_init, "n") <- net_size

nw <- network(el_init, type="edgelist", directed = FALSE, bip=bip_size)

# larger than necessary for bip network but it doesn't hurt anything
wts <- matrix(runif(net_size*net_size), net_size, net_size)

eattr <- wts[el_init]
vattr <- c(rep(seq_len(b1attr_levels), length.out=bip_size), rep(seq_len(b2attr_levels), length.out=net_size - bip_size))
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(seq_len(b1attr_levels*b2attr_levels), b1attr_levels, b2attr_levels)

b1levs <- c(1,2,4,5)
b2levs <- -c(3,6)
levs2 <- c(34,23,11,6,3,9,19,25,28,30)
levs2withblevs <- c(16,4,12,3,19,17,8)

el_off <- el[toggle_off,,drop=FALSE]
el_on <- el[indices[-init_indices],,drop=FALSE]
changes_godfather <- rbind(cbind(el_off, 0), cbind(el_on, wts[el_on]))
changes_godfather <- changes_godfather[sample(seq_len(NROW(changes_godfather))),,drop=FALSE]

el_final <- el[final_indices,,drop=FALSE]
el_final <- cbind(el_final, wts[el_final])
el_final <- el_final[sample(seq_len(NROW(el_final))),,drop=FALSE]

stats <- rep(0, b1attr_levels*b2attr_levels)
types <- m_ind[matrix(vattr[el_final], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + wts[el_final[i,1:2,drop=FALSE]]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="sum", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="sum"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

stats <- rep(0, b1attr_levels*b2attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="nonzero", levels2=TRUE),response="eattr",changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="nonzero"),response="eattr",changes=changes_godfather), ignore_attr = TRUE)

changes_godfather <- changes_godfather[,1:2,drop=FALSE]
expect_equal(stats, ergm.godfather(nw ~ nodemix("vattr", levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[levs2], ergm.godfather(nw ~ nodemix("vattr", levels2 = levs2),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2=TRUE),changes=changes_godfather), ignore_attr = TRUE)
expect_equal(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]], ergm.godfather(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs),changes=changes_godfather), ignore_attr = TRUE)

})


data(faux.mesa.high)
fmh <- faux.mesa.high
set.seed(7)
set.edge.attribute(fmh, "GradeMet", rbinom(203, 6, .5))

fmh.mm.Race <-
  c(mix.Race.Black.Black = 0, mix.Race.Black.Hisp = 8,
    mix.Race.Hisp.Hisp = 53, mix.Race.Black.NatAm = 13,
    mix.Race.Hisp.NatAm = 41, mix.Race.NatAm.NatAm = 46,
    mix.Race.Black.Other = 0, mix.Race.Hisp.Other = 1,
    mix.Race.NatAm.Other = 0, mix.Race.Other.Other = 0,
    mix.Race.Black.White = 5, mix.Race.Hisp.White = 22,
    mix.Race.NatAm.White = 10, mix.Race.Other.White = 0,
    mix.Race.White.White = 4)


test_that("Undirected nodemix() summary with level2 filter by logical matrix", {
  M <- matrix(TRUE, 5, 5)
  M[1,1] <- M[3,2] <- M[2,3] <- FALSE
  s.ab2 <- summary(fmh ~ nodemix("Race", levels2=M))
  expect_equal(s.ab2, fmh.mm.Race[M[upper.tri(M, TRUE)]], ignore_attr=TRUE)
})

test_that("Undirected nodemix() summary with level2 filter by numeric matrix", {
  M <- cbind(2,3)
  s.ab2 <- summary(fmh ~ nodemix("Race", levels2=M))
  expect_equal(s.ab2, fmh.mm.Race[5], ignore_attr=TRUE)
})


data(sampson)
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")

samp.mm.Trinity <- c(
  Trinity.F.F = 8, Trinity.H.F = 12, Trinity.S.F = 12,
  Trinity.F.H = 8, Trinity.H.H = 5, Trinity.S.H = 9,
  Trinity.F.S = 13, Trinity.H.S = 12, Trinity.S.S = 9)

test_that("Directed nodemix() summary with level2 matrix filter", {
  M <- matrix(FALSE, 3, 3)
  M[1,2] <- M[1,3] <- TRUE
  s.ab2 <- summary(samplike ~ nodemix("Trinity", levels2=M))
  expect_equal(s.ab2, samp.mm.Trinity[c(M)], ignore_attr=TRUE)
})

test_that("Undirected nodemix() with level2 character matrix", {
  M <- matrix(NA, 6, 6)
  M[] <- letters[c(abs(row(M) - col(M)))+1]
  M[lower.tri(M,TRUE)] <- NA
  M[which(M == 'b')] <- 'g'

  s <- summary(faux.mesa.high~nodemix("Grade", levels2=M))
  expect_equal(s, c(mix.Grade.c=15, mix.Grade.d=7, mix.Grade.e=2, mix.Grade.f=1, mix.Grade.g=15))
})

test_that("Directed nodemix() with level2 character matrix", {
  M <- matrix(NA, 3, 3)
  M[] <- c(abs(row(M) - col(M)))+1
  M[lower.tri(M, T)] <- M[lower.tri(M, T)]+3
  M[] <- letters[M]

  s <- summary(samplike~nodemix("Trinity", levels2=M))
  expect_equal(s, c(mix.Trinity.b=20, mix.Trinity.c=13, mix.Trinity.d=22, mix.Trinity.e=21, mix.Trinity.f=12))
})

test_that("Bipartite nodemix() with levels2 character matrix", {
  M <- matrix(c('d', 'b', 'b', 'c'), 2, 2)

  s <- summary(mynw~nodemix("names", levels2=M))
  expect_equal(s, c(mix.names.b=13, mix.names.c=1, mix.names.d=9))
})

test_that("Undirected nodemix() with levels2 character matrix with both grouped and ungrouped elements", {
  M <- matrix(NA, 6, 6)
  M[] <- letters[c(abs(row(M)-col(M)))+1]
  M[lower.tri(M,TRUE)] <- NA
  M[1,2] <- M[2,3] <- M[3,4] <- M[1,6] <- ""

  s <- summary(faux.mesa.high~nodemix("Grade", levels2=M))

  ## Expected output constructed as follows;
  M[1,2] <- M[2,3] <- M[3,4] <- M[1,6] <- NA
  grouped <- summary(faux.mesa.high~nodemix("Grade", levels2=M))
  M <- matrix(FALSE, 6,6)
  M[1,2] <- M[2,3] <- M[3,4] <- M[1,6] <- TRUE
  ungrouped <- summary(faux.mesa.high~nodemix("Grade", levels2=M))
  ref <- c(grouped, ungrouped)

  expect_equal(s, ref)
})

test_that("Directed nodemix() with levels2 character matrix with both grouped and ungrouped elements", {
  M <- matrix(c(
    'a', 'b', 'c',
    'd', '', 'b',
    'c', 'd', ''), nrow=3, byrow=TRUE)
  s <- summary(samplike~nodemix("Trinity", levels2=M))

  M[which(M == '')] <- NA
  grouped <- summary(samplike~nodemix("Trinity", levels2=M))
  ungrouped <- summary(samplike~nodemix("Trinity", levels2=is.na(M)))
  expect_equal(s, c(grouped, ungrouped))
})

test_that("Bipartite nodemix() with levels2 character matrix with both grouped and ungrouped elements", {
  M <- matrix(c('', 'b', 'b', ''), 2, 2)

  s <- summary(mynw~nodemix("names", levels2=M))

  M[which(M == '')] <- NA
  grouped <- summary(mynw~nodemix("names", levels2=M))
  ungrouped <- summary(mynw~nodemix("names", levels2=is.na(M)))
  expect_equal(s, c(grouped, ungrouped))
})
