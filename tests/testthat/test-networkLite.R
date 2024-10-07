#  File tests/testthat/test-networkLite.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

## tests are run conditionally on the availability of the networkLite package
if(require("networkLite")) {

  test_that("network and networkLite simulate and summarize formulas equally in ergm", {

    net_size <- 100
    bip_size <- 40

    ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) + gwesp +
              gwnsp(0.3, fixed=TRUE)
    ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp +
                gwnsp(0.3, fixed=TRUE)

    for(directed in list(FALSE, TRUE)) {
      for(bipartite in list(FALSE, bip_size)) {
        if(directed && bipartite) {
          next
        }

        set.seed(0)
        nw <- network.initialize(net_size, directed = directed,
                                 bipartite = bipartite)
        nw %v% "a" <- rep(letters[1:5], length.out = net_size)
        nw %v% "b" <- runif(net_size)

        nwL <- as.networkLite(nw)

        coef <- c(-4, 1, 1.5, 0.5, -1, 0.5)

        set.seed(0)
        nw_1 <- simulate(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                         coef = coef, output = "network")
        set.seed(0)
        nwL_1 <- simulate(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                          coef = coef, output = "network")
        expect_s3_class(nwL_1, "networkLite")

        expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
        if(directed) {
          expect_identical(summary(ffdir, basis = nw_1),
                           summary(ffdir, basis = nwL_1))
        } else {
          expect_identical(summary(ffundir, basis = nw_1),
                           summary(ffundir, basis = nwL_1))
        }

        set.seed(0)
        nw_2 <- simulate(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                         coef = coef, output = "network")
        set.seed(0)
        nwL_2 <- simulate(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                          coef = coef, output = "network")
        expect_s3_class(nwL_2, "networkLite")

        expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
        if(directed) {
          expect_identical(summary(ffdir, basis = nw_2),
                           summary(ffdir, basis = nwL_2))
        } else {
          expect_identical(summary(ffundir, basis = nw_2),
                           summary(ffundir, basis = nwL_2))
        }
      }
    }
  })

  test_that("network and networkLite simulate equally in san", {

    net_size <- 100
    bip_size <- 40

    ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) + gwesp +
              gwnsp(0.3, fixed=TRUE)
    ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp +
                gwnsp(0.3, fixed=TRUE)

    for(directed in list(FALSE, TRUE)) {
      for(bipartite in list(FALSE, bip_size)) {
        if(directed && bipartite) {
          next
        }

        set.seed(0)
        nw <- network.initialize(net_size, directed = directed,
                                 bipartite = bipartite)
        nw %v% "a" <- rep(letters[1:5], length.out = net_size)
        nw %v% "b" <- runif(net_size)

        nwL <- as.networkLite(nw)

        set.seed(0)
        nw_1 <- san(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                    target.stats = c(1000, 500, 300, 200, 600, 1500))
        set.seed(0)
        nwL_1 <- san(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                     target.stats = c(1000, 500, 300, 200, 600, 1500))
        expect_s3_class(nwL_1, "networkLite")

        expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
        if(directed) {
          expect_identical(summary(ffdir, basis = nw_1),
                           summary(ffdir, basis = nwL_1))
        } else {
          expect_identical(summary(ffundir, basis = nw_1),
                           summary(ffundir, basis = nwL_1))
        }

        set.seed(0)
        nw_2 <- san(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                    target.stats = c(800, 400, 200, 100, 600, 1200))
        set.seed(0)
        nwL_2 <- san(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b),
                     target.stats = c(800, 400, 200, 100, 600, 1200))
        expect_s3_class(nwL_2, "networkLite")

        expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
        if(directed) {
          expect_identical(summary(ffdir, basis = nw_2),
                           summary(ffdir, basis = nwL_2))
        } else {
          expect_identical(summary(ffundir, basis = nw_2),
                           summary(ffundir, basis = nwL_2))
        }
      }
    }
  })

  test_that("network and networkLite fit and simulate equal missing-data ergms", {

    net_size <- 50
    bip_size <- 20

    for(directed in list(FALSE, TRUE)) {
      for(bipartite in list(FALSE, bip_size)) {
        if(directed && bipartite) {
          next
        }
        if(directed) {
          ergm_formula <- ~edges + odegree(1) + absdiff("age")
        } else {
          ergm_formula <- ~edges + degree(1) + absdiff("age")
        }
        set.seed(0)
        nwL <- networkLite(net_size, directed = directed, bipartite = bipartite)
        nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL)/10)
        nwL %v% "age" <- runif(net_size)
        na <- sample(c(FALSE,TRUE),network.edgecount(nwL),TRUE)

        set.seed(0)
        eL <- ergm(ergm_formula, basis = nwL,
                   control = list(MCMLE.effectiveSize = NULL))
        set.edge.attribute(nwL, "na", na)
        set.seed(0)
        eLna <- ergm(ergm_formula, basis = nwL,
                     control = list(MCMLE.effectiveSize = NULL))
        eL2 <- simulate(eLna)
        expect_s3_class(eL2, "networkLite")

        set.seed(0)
        nw <- network.initialize(net_size, directed = directed,
                                 bipartite = bipartite)
        nw <- san(nw ~ edges, target.stats = network.dyadcount(nw)/10)
        nw %v% "age" <- runif(net_size)
        na <- sample(c(FALSE,TRUE),network.edgecount(nw),TRUE)

        set.seed(0)
        e <- ergm(ergm_formula, basis = nw,
                  control = list(MCMLE.effectiveSize = NULL))
        set.edge.attribute(nw, "na", na)
        set.seed(0)
        ena <- ergm(ergm_formula, basis = nw,
                    control = list(MCMLE.effectiveSize = NULL))
        e2 <- simulate(ena)

        expect_equal(coef(e), coef(eL))
        expect_equal(coef(ena), coef(eLna))
        expect_equal(as.edgelist(e2), as.edgelist(eL2))
        expect_equal(as.edgelist(e2, attrname = "na"),
                     as.edgelist(eL2, attrname = "na"))
      }
    }
  })

  test_that("network and networkLite fit and simulate equal valued ergms", {

    net_size <- 50
    bip_size <- 20

    for(directed in list(FALSE, TRUE)) {
      for(bipartite in list(FALSE, bip_size)) {
        if(directed && bipartite) {
          next
        }

        set.seed(0)
        nwL <- networkLite(net_size, directed = directed,
                           bipartite = bipartite)
        nwL <- san(nwL ~ edges, target.stats = network.dyadcount(nwL))
        nwL %v% "age" <- runif(net_size)
        set.edge.attribute(nwL, "w", runif(network.edgecount(nwL)))
        eL <- ergm(nwL ~ absdiff("age"), response = "w", reference = ~Unif(0,1),
                   control = list(MCMLE.effectiveSize = NULL))
        eL2 <- simulate(eL)
        expect_s3_class(eL2, "networkLite")

        set.seed(0)
        nw <- network.initialize(net_size, directed = directed,
                                 bipartite = bipartite)
        nw <- san(nw ~ edges, target.stats = network.dyadcount(nw))
        nw %v% "age" <- runif(net_size)
        set.edge.attribute(nw, "w", runif(network.edgecount(nw)))
        e <- ergm(nw ~ absdiff("age"), response = "w", reference = ~Unif(0,1),
                  control = list(MCMLE.effectiveSize = NULL))
        e2 <- simulate(e)

        expect_equal(coef(e), coef(eL))
        expect_equal(as.edgelist(e2, attrname = "w"),
                     as.edgelist(eL2, attrname = "w"))
      }
    }
  })

}
