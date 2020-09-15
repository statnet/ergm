#  File tests/nodemix.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
library(ergm)
# Test undirected network
data(faux.mesa.high)
m <- matrix(c(75, 0, 0, 1, 1, 1, 0, 33, 2, 4, 2, 1,
 0, 2, 23, 7, 6, 4, 1, 4, 7, 9, 1, 5, 1, 2, 6, 1, 17, 5,
 1, 1, 4, 5, 5, 6), 6, 6)  # Correct answer!
faux.mesa.high %e% "eattr" <- rep(1, network.edgecount(faux.mesa.high))
if (max(abs(m[upper.tri(m, diag=T)] - summary(faux.mesa.high ~ nodemix("Grade")))) > 0)
  stop ("binary nodemix failed test on undirected network faux.mesa.high")
if (max(abs(m[upper.tri(m, diag=T)] - summary(faux.mesa.high ~ nodemix("Grade", form="nonzero"), response="eattr"))) > 0)
  stop ("weighted nodemix failed test on undirected network faux.mesa.high")


# directed network
data(sampson)
grpord<-c("Turks","Loyal","Outcasts")
m2 <- matrix(c(30, 9, 7, 5, 23, 1, 1, 2, 10), 3, 3,dimnames=list(From=grpord,To=grpord)) # Correct answer!
mixnames<-t(sapply(grpord,function(from) sapply(grpord,function(to) paste("mix.group",from,to,sep="."))))
samplike %e% "eattr" <- rep(1, network.edgecount(samplike))
if (!all(c(m2) == summary(samplike ~ nodemix("group"))[c(mixnames)]))
  stop ("binary nodemix failed test on directed network samplike")
if (!all(c(m2) == summary(samplike ~ nodemix("group", form="nonzero"), response="eattr")[c(mixnames)]))
  stop ("weighted nodemix failed test on directed network samplike")

# bipartite network
el <- cbind( c(17, 20, 22, 26, 19, 24, 16, 22, 18, 23, 28, 20,
               22, 23, 17, 21, 25, 21, 27, 16, 19, 18, 23),
           c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 10, 
             10, 11, 11))
mynw <- network(el, bipartite=15, directed=FALSE) 
mynw %v% "names" <- rep(letters[1:3], c(10,10,8))
m3 <- matrix(c(9, 1, 12, 1), 2, 2) # Correct answer!
mynw %e% "eattr" <- rep(1, network.edgecount(mynw))
if (max(abs(as.vector(m3) - summary(mynw ~ nodemix("names")))) > 0)
  stop ("binary nodemix failed test on bipartite network")
if (max(abs(as.vector(m3) - summary(mynw ~ nodemix("names", form="nonzero"), response="eattr"))) > 0)
  stop ("weighted nodemix failed test on bipartite network")

net_size <- 1000
bip_size <- 300
attr_levels <- 6
b1attr_levels <- 5
b2attr_levels <- 7

## undirected ##

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
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr")))

stats <- rep(0, attr_levels*(attr_levels + 1)/2)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs))))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs))))


## directed ##

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
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="sum"),response="eattr")))

stats <- rep(0, attr_levels*attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs, form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))] == summary(nw ~ nodemix("vattr", levels = levs))))
stopifnot(all(stats[sort(unique(c(m_ind[levs,levs])))[levs2withlevs]] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = levs2withlevs))))

## bipartite ##

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
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="sum"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="sum"),response="eattr")))

stats <- rep(0, b1attr_levels*b2attr_levels)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs, form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs))))
stopifnot(all(stats[sort(unique(c(m_ind[b1levs,b2levs])))[levs2withblevs]] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = levs2withblevs))))
