#  File tests/testthat/test-metrics.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
attach(MLE.tools)

theta0err<- 1 # Perturbation in the initial values

n<-20 # Number of nodes

d<-.1 # Density
ms <-c(0,.05) # Missingness rates
metrics <- c("naive", "lognormal", "median", "logtaylor", "EF.Likelihood")
EF_warned <- FALSE

run.metric.test<-function(y){
  truth<-edges.theta(y)

  for(metric in metrics){
    test_that(paste0("Metric test for: ", format(metric), ", n = ", n, ", naive density = ", format(network.edgecount(y)/network.dyadcount(y)), ", missing fraction = ", format(network.naedgecount(y)/network.dyadcount(y)), "."), {
      torun <- quote(mcmcfit <- ergm(y~edges, control=control.ergm(force.main=TRUE, init=truth+theta0err, MCMLE.metric=metric),eval.loglik=FALSE, verbose=FALSE))
      if (network.naedgecount(y) && metric == "logtaylor")
        expect_error(eval(torun), "Metric 'logtaylor' is not implemented.*")
      else {
        # This is not 100% safe on the off chance that this warning
        # gets invoked elsewhere in the tests.
        if (metric == "EF.Likelihood" && !EF_warned) {
          expect_warning(eval(torun), "Metric 'EF.Likelihood' has been deprecated in favor of 'naive'.*")
          EF_warned <<- TRUE
        } else {
          eval(torun)
        }
        expect_within_mc_err(mcmcfit, truth)
      }
    })
  }
}

for(m in ms){
  set.seed(123)
  y<-mk.missnet(n, d, m, TRUE, FALSE)
  run.metric.test(y)
}

detach(MLE.tools)
