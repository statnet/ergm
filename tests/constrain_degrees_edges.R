#  File tests/constrain_degrees_edges.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# Also, while we are at it, test the options setting in .onLoad().
options(ergm.eval.loglik=FALSE) # .onLoad() should not clobber this.

library(ergm)

# Check that options are either set to default or preserved.
stopifnot(getOption("ergm.eval.loglik")==FALSE)
stopifnot(getOption("ergm.loglik.warn_dyads")==TRUE)

nsim <- 100
n <- 50
m <- 10
d <- 0.1

od <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 1, sum)
id <- function(nw) apply(as.matrix(nw, matrix.type="adjacency"), 2, sum)
e <- function(nw) network.edgecount(nw)

###### Directed
y0 <- as.network(n, density=d, directed=TRUE)

### Outdegrees
ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees, coef=rep(0,n*2), nsim=nsim, output="stats")
stopifnot(all(sweep(ys[,1:n], 2, od(y0))==0), any(sweep(ys[,-(1:n)], 2, id(y0))!=0)) # Outdegrees shouldn't vary but indegrees should.

### Indegrees
ys <- simulate(y0~receiver(nodes=TRUE)+sender(nodes=TRUE), constraints=~idegrees, coef=rep(0,n*2), nsim=nsim, output="stats")
stopifnot(all(sweep(ys[,1:n], 2, id(y0))==0), any(sweep(ys[,-(1:n)], 2, od(y0))!=0)) # Indegrees shouldn't vary but outdegrees should.

### Both in- and outdegrees
ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~degrees, coef=rep(0,n*2), nsim=nsim, output="stats")
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~odegrees+idegrees, coef=rep(0,n*2), nsim=nsim, output="stats")
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

### Edges
ys <- simulate(y0~sender(nodes=TRUE)+receiver(nodes=TRUE), constraints=~edges, coef=rep(0,n*2), nsim=nsim, output="stats")
stopifnot(all(e(y0)==rowSums(ys[,1:n])), all(e(y0)==rowSums(ys[,-(1:n)])),
          any(sweep(ys[,1:n], 2, od(y0))!=0), any(sweep(ys[,-(1:n)], 2, id(y0))!=0)) # Edges shouldn't vary, but in- and out-degrees should.

###### Undirected
y0 <- as.network(n, density=d, directed=FALSE)

### Degrees
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~degrees, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(sweep(ys, 2, od(y0))==0))

### Edges
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~edges, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(e(y0)==rowSums(ys)/2),
          any(sweep(ys, 2, od(y0))!=0)) # Edges shouldn't vary, but degrees should.

###### Bipartite undirected
y0 <- as.network(n, density=d, directed=FALSE, bipartite=m)

### B1degrees
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b1degrees, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(sweep(ys[,1:m], 2, od(y0))==0), any(sweep(ys[,-(1:m)], 2, id(y0))!=0))

### B2degrees
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b2degrees, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(sweep(ys[,-(1:m)], 2, id(y0))==0), any(sweep(ys[,1:m], 2, od(y0))!=0))

### Both B1 and B2 degrees
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~degrees, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

ys <- simulate(y0~sociality(nodes=TRUE), constraints=~b1degrees+b2degrees, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(sweep(ys, 2, c(od(y0),id(y0)))==0))

### Edges
ys <- simulate(y0~sociality(nodes=TRUE), constraints=~edges, coef=rep(0,n), nsim=nsim, output="stats")
stopifnot(all(e(y0)==rowSums(ys)/2),
          any(sweep(ys, 2, c(od(y0),id(y0)))!=0)) # Edges shouldn't vary, but degrees should.
