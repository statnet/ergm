#  File tests/testthat/test-constraints.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

set.seed(0)

net1 <- network.initialize(10,directed=FALSE)
net1[,] <- 1
absent <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]
net1[absent] <- 0
present <- as.edgelist(net1)[sample.int(network.edgecount(net1), 2), ]
fixed <- rbind(present, absent)

net1[as.edgelist(net1)[sample.int(network.edgecount(net1), round(network.edgecount(net1)/2)), ]] <- 0
net1[present] <- 1

test_that("fixedas(present, absent)", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent))
  s1 <- simulate(t1, 100)

  # check if all the simulated network have 'present' edges
  expect_true(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

  # check if all the simulated network do not have 'absent' edges
  expect_true(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))
})


test_that("fixedas(fixed.dyads)", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(fixed))
  s1 <- simulate(t1, 100)

  # check that fixed edges are identical between simulated networks and the original network
  expect_true(all(sapply(s1,function(x) identical(x[fixed], net1[fixed]))))
})


test_that("only present", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(present = present))
  s1 <- simulate(t1,100)
  expect_true(all(sapply(s1,function(x)as.data.frame(t(present)) %in% as.data.frame(t(as.edgelist(x))))))

  # Also test for inconsistent constraint.
  expect_error(ergm(net1~edges, constraint = ~fixedas(present = absent)),
               "In constraint 'fixedas' in package 'ergm': Edges constrained to be present are absent in the LHS network.")
})

test_that("only absent", {
  t1 <- ergm(net1~edges, constraint = ~fixedas(absent = absent))
  s1 <- simulate(t1, 100)
  expect_true(all(!sapply(s1,function(x)as.data.frame(t(absent)) %in% as.data.frame(t(as.edgelist(x))))))

  # Also test for inconsistent constraint.
  expect_error(ergm(net1~edges, constraint = ~fixedas(absent = present)),
               "In constraint 'fixedas' in package 'ergm': Edges constrained to be absent are present in the LHS network.")
})

present <- as.network(present, matrix.type = "edgelist", directed = FALSE)
absent <- as.network(absent, matrix.type = "edgelist", directed = FALSE)

test_that("fixedas with network input", {
  expect_warning(t1 <- ergm(net1~edges, constraint = ~fixedas(present = present, absent = absent)),
                 "^In constraint 'fixedas' in package 'ergm': Network size of argument\\(s\\) 'present' and 'absent' differs from that of the response network\\..*")
  expect_warning(s1 <- simulate(t1, 100),
                 "^In constraint 'fixedas' in package 'ergm': Network size of argument\\(s\\) 'present' and 'absent' differs from that of the response network\\..*")

  expect_true(all(sapply(s1,function(x)as.data.frame(t(as.edgelist(present))) %in% as.data.frame(t(as.edgelist(x))))))
  expect_true(all(!sapply(s1,function(x)as.data.frame(t(as.edgelist(absent))) %in% as.data.frame(t(as.edgelist(x))))))
})

net1 <- network(10,directed=FALSE,density=0.5)
fdel <- matrix(sample(2:9,8,replace=FALSE),4,2)

for(free.dyads in list(
                    fdel,
                    fdnw <- as.network(structure(fdel, n = 10), directed = FALSE),
                    fd <- as.rlebdm(fdnw)
                  )){
  test_that(sprintf("fixallbut with %s input", class(free.dyads)[1]), {
    t1 <- ergm(net1~edges, constraint = ~fixallbut(free.dyads = free.dyads))
    s1 <- simulate(t1, 100)

    fixed.dyads <- as.edgelist(!update(net1,fdel,matrix.type="edgelist"))
    fixed.dyads.state <- net1[fixed.dyads]

    expect_true(all(sapply(s1,function(x) all.equal(x[fixed.dyads],fixed.dyads.state))))
  })
}


test_that("constraint conflict is detected", {
  data(florentine)
  conwarn <- "^The specified model's sample space constraint holds statistic\\(s\\) edges  constant. They will be ignored.$"
  dyadwarn <- "^The number of observed dyads in this network is ill-defined due to complex constraints on the sample space..*$"
  
  ergm(flomarriage~edges, constraints = ~edges) |>
    expect_warning(conwarn) |>
    expect_warning(dyadwarn) |>
    expect_warning(dyadwarn)

  (fit <- ergm(flomarriage~edges + triangle, constraints = ~degrees)) |>
    expect_warning(conwarn) |>
    expect_warning(dyadwarn) |>
    expect_warning(dyadwarn)

  expect_equal(coef(fit)[1],0, ignore_attr=TRUE)
})

test_that("ChangeStats() constraining a dyad-dependent statistic to 0", {
  S <- 100
  data(florentine)

  flo0 <- flomarriage
  flo0[,] <- 0

  sim <- simulate(flo0 ~ edges + triangle, coef = c(0, 1), constraints = ~ChangeStats(~triangle), output = "stats", nsim = S)
  expect_equal(sim[, 2], numeric(S))
})

test_that("ChangeStats() constraining a dyad-independent statistic falls back to Dyads()", {
  S <- 100
  data(sampson)

  expect_message(sim <- simulate(samplike ~ edges, coef = 0, monitor = ~nodemix("group", levels2 = TRUE), constraints = ~ChangeStats(~nodematch("group")), output = "stats", nsim = S),
                 ".*constraint formula is dyad-independent; falling back.*")

  expect_equal(apply(sim, 2, sd)[c(2,6,10)], numeric(3), ignore_attr = TRUE)
  expect_true(all(apply(sim, 2, sd)[-c(2,6,10)] > 0))
})

test_that("ChangeStats() constraining a dyad-independent statistic enforces", {
  S <- 100
  data(sampson)

  sim <- expect_silent(simulate(samplike ~ edges, coef = 0, monitor = ~nodemix("group", levels2 = TRUE), constraints = ~ChangeStats(~nodematch("group"), FALSE), output = "stats", nsim = S))

  expect_equal(apply(sim, 2, sd)[c(2,6,10)], numeric(3), ignore_attr = TRUE)
  expect_true(all(apply(sim, 2, sd)[-c(2,6,10)] > 0))
})
