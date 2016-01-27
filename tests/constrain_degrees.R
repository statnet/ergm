#  File tests/constrain_degrees.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
library(ergm)

nsim <- 100
n <- 50
m <- 10
d <- 0.1

od <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 1, sum)
id <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 2, sum)

###### Directed
y0 <- as.network(n, density=d, directed=TRUE)

### Outdegrees
ys <- simulate(y0~sender(base=0)+receiver(base=0), constraints=~odegrees, coef=rep(0,n*2), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys[,1:n], 2, od(y0))==0), any(sweep(ys[,-(1:n)], 2, id(y0))!=0)) # Outdegrees shouldn't vary but indegrees should.

### Indegrees
ys <- simulate(y0~receiver(base=0)+sender(base=0), constraints=~idegrees, coef=rep(0,n*2), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys[,1:n], 2, id(y0))==0), any(sweep(ys[,-(1:n)], 2, od(y0))!=0)) # Indegrees shouldn't vary but outdegrees should.

### Both in- and outdegrees
ys <- simulate(y0~sender(base=0)+receiver(base=0), constraints=~degrees, coef=rep(0,n*2), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

ys <- simulate(y0~sender(base=0)+receiver(base=0), constraints=~odegrees+idegrees, coef=rep(0,n*2), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

###### Undirected
y0 <- as.network(n, density=d, directed=FALSE)

### Degrees
ys <- simulate(y0~sociality(base=0), constraints=~degrees, coef=rep(0,n), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys, 2, od(y0))==0))

###### Bipartite undirected
y0 <- as.network(n-m, density=d, directed=FALSE, bipartite=m)

### B1degrees
ys <- simulate(y0~sociality(base=0), constraints=~b1degrees, coef=rep(0,n), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys[,1:m], 2, od(y0))==0), any(sweep(ys[,-(1:m)], 2, id(y0))!=0))

### B2degrees
ys <- simulate(y0~sociality(base=0), constraints=~b2degrees, coef=rep(0,n), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys[,-(1:m)], 2, id(y0))==0), any(sweep(ys[,1:m], 2, od(y0))!=0))

### Both B1 and B2 degrees
ys <- simulate(y0~sociality(base=0), constraints=~degrees, coef=rep(0,n), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

ys <- simulate(y0~sociality(base=0), constraints=~b1degrees+b2degrees, coef=rep(0,n), nsim=nsim, statsonly=TRUE)
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))
