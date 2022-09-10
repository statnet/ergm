#  File tests/testthat/test-parallel.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

library(statnet.common)
opttest({

data(florentine)

for(type in c("PSOCK")){
  test_that(paste0("parallel ", type), {
    gest <- ergm(flomarriage ~ edges + absdiff("wealth"),
                 eval.loglik=TRUE,
                 control=control.ergm(MCMC.burnin=1000, MCMC.interval=10, MCMLE.maxit=2, MCMC.samplesize=1000, force.main=TRUE,
                   parallel=2, parallel.type=type))

    print(summary(gest))
    mcmc.diagnostics(gest)

    # FIXME: Set seeds and replace with actual values?
    sim.STAT.SEQ <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), output="stats", sequential=TRUE)
    expect_equal(nrow(sim.STAT.SEQ), 5)

    sim.STAT.seq <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), output="stats", sequential=FALSE)
    expect_equal(nrow(sim.STAT.seq), 5)

    sim.stat.SEQ <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), output="network", sequential=TRUE)
    expect_equal(length(sim.stat.SEQ), 5)

    sim.stat.seq <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), output="network", sequential=FALSE)
    expect_equal(length(sim.stat.seq), 5)

    if(exists("cl")){
      stopCluster(cl)
      rm(cl)
    }
  })
}

# pre-made SOCK cluster
test_that("pre-made SOCK cluster", {
  library(parallel)
  clus <- makeCluster(2, type='PSOCK')
  clus
  data(florentine)
  t0 <- proc.time()
  fauxmodel.01 <- ergm(flomarriage ~ edges + isolates + gwesp(0.2, fixed=T),
                       control=control.ergm(parallel=clus,
                                            MCMLE.maxit=100))
  proc.time() - t0

  sim.test <- simulate(network(1000) ~ edges +
                             gwesp(0.5, fixed = TRUE) + mutual + istar(2),
                           coef = c(-5, 1.3, 1.5, -0.5), nsim = 4,
                           control = control.simulate.formula(parallel = clus,
                                                              parallel.type = "PSOCK"))

  stopCluster(clus)
})

}, "parallel")


if(inherits(try(get.MT_terms(), silent=TRUE),"try-error")){
  message("Skipping OpenMP test. This package installation was built without OpenMP support.")
}else{
  test_that("OpenMP", {
    library(ergm)
    data(florentine)
    set.seed(0)
    sim.ser <- simulate(flomarriage~edges+triangle, nsim=100, control=control.simulate(MCMC.burnin=1, MCMC.interval=1), output="stats")

    prev <- set.MT_terms(2)
    set.seed(0)
    sim.par <- simulate(flomarriage~edges+triangle, nsim=100, control=control.simulate(MCMC.burnin=1, MCMC.interval=1), output="stats")
    set.MT_terms(prev)

    expect_true(all.equal(sim.ser,sim.par))
  })
}
