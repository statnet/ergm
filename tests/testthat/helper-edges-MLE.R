#  File tests/testthat/helper-edges-MLE.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
logit<-function(p) log(p/(1-p))

MLE.tools <- new.env()

MLE.tools$mk.missnet<-function(n,d,m,directed=TRUE,bipartite=FALSE){
  y<-network.initialize(n, directed=directed, bipartite=bipartite)
  y<-simulate(y~edges, coef=logit(d), control=control.simulate(MCMC.burnin=2*n^2))
  if(m>0){
    y.miss<-simulate(y~edges, coef=logit(m))
    y[as.edgelist(y.miss)]<-NA
  }
  y
}

MLE.tools$edges.theta<-function(y){
  e<-network.edgecount(y)
  d<-network.dyadcount(y)
  logit(e/d)
}

MLE.tools$edges.llk<-function(y, theta=NULL, e=NULL){
  e<-NVL(e, network.edgecount(y))
  d<-network.dyadcount(y)
  NVL(theta) <- logit(e/d)
  e*theta - d*log1p(exp(theta))
}
