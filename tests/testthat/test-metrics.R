#  File tests/testthat/test-metrics.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
theta0err<- 1 # Perturbation in the initial values
tolerance<-5 # Result must be within 5*MCMCSE of truth.

n<-20 # Number of nodes

d<-.1 # Density
ms <-c(0,.05) # Missingness rates
metrics <- c("naive", "lognormal", "Median.Likelihood")

logit<-function(p) log(p/(1-p))

mk.missnet<-function(n,d,m,directed=TRUE,bipartite=FALSE){
  y<-network.initialize(n, directed=directed, bipartite=bipartite)
  y<-simulate(y~edges, coef=logit(d), control=control.simulate(MCMC.burnin=2*n^2))
  if(m>0){
    y.miss<-simulate(y~edges, coef=logit(m))
    y[as.edgelist(y.miss)]<-NA
  }
  y
}

correct.edges.theta<-function(y){
  e<-network.edgecount(y)
  d<-network.dyadcount(y)

  logit(e/d)
}

run.metric.test<-function(y){
  truth<-correct.edges.theta(y)

  for(metric in metrics){
    test_that(paste0("Metric test for: ", metric, ", n = ", n, ", naive density = ", network.edgecount(y)/network.dyadcount(y), ", missing fraction = ", network.naedgecount(y)/network.dyadcount(y), "."), {
      mcmcfit<-ergm(y~edges, control=control.ergm(force.main=TRUE, init=truth+theta0err, MCMLE.metric=metric),eval.loglik=FALSE, verbose=FALSE)
      mcmcOK<-abs(truth-coef(mcmcfit))/sqrt(diag(vcov(mcmcfit, source="estimation")))
      expect_lt(mcmcOK, tolerance)
    })
  }
}

for(m in ms){
  set.seed(123)
  y<-mk.missnet(n, d, m, TRUE, FALSE)
  run.metric.test(y)
}
