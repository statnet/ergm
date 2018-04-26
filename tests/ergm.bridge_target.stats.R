#  File tests/ergm.bridge_target.stats.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
set.seed(1)
logit <- function(p) log(p/(1-p))
expit <- function(e) exp(e)/(1+exp(e))
l <- function(nw, ets=NULL, theta=NULL){
  e <- NVL(ets, network.edgecount(nw))
  d <- network.dyadcount(nw)
  theta <- NVL(theta, logit(e/d))

  e*log(expit(theta)) + (d-e)*log(expit(-theta))
}

data(florentine)
y <- flomarriage

# Attainable t.s..
ts <- 3
print(llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts))))
print(llk <- l(y,ts))
stopifnot(isTRUE(all.equal(llk,llk.ergm, check.attributes=FALSE)))

# Unattainable t.s..
ts <- 3.5
print(llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts))))
print(llk <- l(y,ts))
stopifnot(isTRUE(all.equal(llk,llk.ergm, check.attributes=FALSE,tolerance=0.05)))
}, "bridge sampling with target stats")
