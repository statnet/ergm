#  File tests/constraint_generalized.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

library(ergm)

# fixedas
net1 <- network.initialize(10,directed=FALSE)
net1[,] <- 1
absent <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]
net1[absent] <- 0
present <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]

net1[as.edgelist(net1)[sample.int(network.edgecount(net1), round(network.edgecount(net1)/2)), ]] <- 0
net1[present] <- 1

t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent))
s1 <- simulate(t1, 100)

# check if all the simulated network have 'present' edges
stopifnot(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

# check if all the simulated network do not have 'absent' edges
stopifnot(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))

# only present
t1 <- ergm(net1~edges, constraint = ~fixedas(present = present))
s1 <- simulate(t1,100)
stopifnot(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

# only absent
t1 <- ergm(net1~edges, constraint = ~fixedas(absent = absent))
s1 <- simulate(t1, 100)
stopifnot(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))

# input is network
present <- as.network(present, matrix.type = "edgelist", directed = FALSE)
absent <- as.network(absent, matrix.type = "edgelist", directed = FALSE)

t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent))
s1 <- simulate(t1, 100)

stopifnot(all(sapply(s1,function(x)as.data.frame(t(as.edgelist(present))) %in% as.data.frame(t(as.edgelist(x))))))
stopifnot(all(!sapply(s1,function(x)as.data.frame(t(as.edgelist(absent))) %in% as.data.frame(t(as.edgelist(x))))))

# fixallbut
net1 <- network(10,directed=FALSE,density=0.5)
free.dyads <- as.edgelist(matrix(sample(1:10,8,replace=F),4,2),n=10,directed=FALSE)

t1 <- ergm(net1~edges, constraint = ~fixallbut(free.dyads = free.dyads))
s1 <- simulate(t1, 100)

fixed.dyads <- as.edgelist(!update(net1,free.dyads,matrix.type="edgelist"))
fixed.dyads.state <- net1[fixed.dyads]

stopifnot(all(sapply(s1,function(x) all.equal(x[fixed.dyads],fixed.dyads.state))))
