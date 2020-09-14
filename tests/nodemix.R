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

## extend below to binary as well, and (b1/b2)levels + levels2 ##

## undirected ##

nw <- network.initialize(100, dir=FALSE)
nw <- san(nw ~ edges, target.stats=100)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- rep(1:3, 100)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(c(1,2,4,2,3,5,4,5,6),3,3)
levs <- -2
levs2 <- c(1,5,3,6)

stats <- rep(0, 6)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[c(1,4,6)] == summary(nw ~ nodemix("vattr", levels = levs, form="sum"),response="eattr")))
stopifnot(all(stats[c(6,1)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(3,1), form="sum"),response="eattr")))

stats <- rep(0, 6)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(1,4,6)] == summary(nw ~ nodemix("vattr", levels = levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(6,1)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(3,1), form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[c(1,4,6)] == summary(nw ~ nodemix("vattr", levels = levs))))
stopifnot(all(stats[c(6,1)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(3,1)))))


## directed ##

nw <- network.initialize(100, dir=TRUE)
nw <- san(nw ~ edges, target.stats=100)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- rep(1:3, 100)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(1:9,3,3)
levs <- 3:2
levs2 <- c(1,5,3,7)

stats <- rep(0, 9)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[c(9,8,6,5)] == summary(nw ~ nodemix("vattr", levels = levs, form="sum"),response="eattr")))
stopifnot(all(stats[c(8,5,9)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(2,4,1), form="sum"),response="eattr")))

stats <- rep(0, 9)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(9,8,6,5)] == summary(nw ~ nodemix("vattr", levels = levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(8,5,9)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(2,4,1), form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[c(9,8,6,5)] == summary(nw ~ nodemix("vattr", levels = levs))))
stopifnot(all(stats[c(8,5,9)] == summary(nw ~ nodemix("vattr", levels = levs, levels2 = c(2,4,1)))))

## bipartite ##

nw <- network.initialize(100, dir=FALSE, bip=30)
nw <- san(nw ~ edges, target.stats=100)
el <- as.edgelist(nw)
eattr <- runif(network.edgecount(nw))
vattr <- rep(1:3, 100)
nw %e% "eattr" <- eattr
nw %v% "vattr" <- vattr
m_ind <- matrix(1:9,3,3)
b1levs <- c(3,1)
b2levs <- -3
levs2 <- c(1,5,3,7)

stats <- rep(0, 9)
types <- m_ind[matrix(vattr[el], ncol=2)]
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + eattr[i]
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="sum"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="sum"),response="eattr")))
stopifnot(all(stats[c(3,1,6,4)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="sum"),response="eattr")))
stopifnot(all(stats[c(4,1,3)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = c(4,2,1), form="sum"),response="eattr")))

stats <- rep(0, 9)
for(i in seq_along(types)) stats[types[i]] <- stats[types[i]] + 1
stopifnot(all(stats == summary(nw ~ nodemix("vattr", form="nonzero"),response="eattr")))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(3,1,6,4)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, form="nonzero"),response="eattr")))
stopifnot(all(stats[c(4,1,3)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = c(4,2,1), form="nonzero"),response="eattr")))

stopifnot(all(stats == summary(nw ~ nodemix("vattr"))))
stopifnot(all(stats[levs2] == summary(nw ~ nodemix("vattr", levels2 = levs2))))
stopifnot(all(stats[c(3,1,6,4)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs))))
stopifnot(all(stats[c(4,1,3)] == summary(nw ~ nodemix("vattr", b1levels = b1levs, b2levels = b2levs, levels2 = c(4,2,1)))))
