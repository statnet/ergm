#  File tests/testthat/test-constraints.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

net1 <- network.initialize(10,directed=FALSE)
net1[,] <- 1
absent <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]
net1[absent] <- 0
present <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]

net1[as.edgelist(net1)[sample.int(network.edgecount(net1), round(network.edgecount(net1)/2)), ]] <- 0
net1[present] <- 1

test_that("fixedas", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent))
  s1 <- simulate(t1, 100)

  # check if all the simulated network have 'present' edges
  expect_true(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

  # check if all the simulated network do not have 'absent' edges
  expect_true(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))
})

test_that("only present", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(present = present))
  s1 <- simulate(t1,100)
  expect_true(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))
})

test_that("only absent", {
t1 <- ergm(net1~edges, constraint = ~fixedas(absent = absent))
s1 <- simulate(t1, 100)
expect_true(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))
})

present <- as.network(present, matrix.type = "edgelist", directed = FALSE)
absent <- as.network(absent, matrix.type = "edgelist", directed = FALSE)

test_that("fixedas with network input", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent))
  s1 <- simulate(t1, 100)

  expect_true(all(sapply(s1,function(x)as.data.frame(t(as.edgelist(present))) %in% as.data.frame(t(as.edgelist(x))))))
  expect_true(all(!sapply(s1,function(x)as.data.frame(t(as.edgelist(absent))) %in% as.data.frame(t(as.edgelist(x))))))
})

test_that("fixallbut with network input", {
  net1 <- network(10,directed=FALSE,density=0.5)
  free.dyads <- matrix(sample(2:9,8,replace=FALSE),4,2)

  t1 <- ergm(net1~edges, constraint = ~fixallbut(free.dyads = free.dyads))
  s1 <- simulate(t1, 100)

  fixed.dyads <- as.edgelist(!update(net1,free.dyads,matrix.type="edgelist"))
  fixed.dyads.state <- net1[fixed.dyads]

  expect_true(all(sapply(s1,function(x) all.equal(x[fixed.dyads],fixed.dyads.state))))
})

test_that("constraint conflict is detected", {
  data(florentine)
  conwarn <- "^The specified model's sample space constraint holds statistic\\(s\\) edges  constant. They will be ignored.$"
  dyadwarn <- "^The number of observed dyads in this network is ill-defined due to complex constraints on the sample space..*$"
  
  expect_warning(ergm(flomarriage~edges, constraints = ~edges), conwarn)

  expect_warning(expect_warning(expect_warning(fit <- ergm(flomarriage~edges + triangle, constraints = ~degrees),
                                               conwarn), dyadwarn), dyadwarn)

  expect_equal(coef(fit)[1],0, ignore_attr=TRUE)
})
