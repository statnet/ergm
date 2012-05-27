#  File ergm/tests/dynamic_MLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)

tolerance<-0.05
n<-10
theta<--1.5

logit<-function(p) log(p/(1-p))

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(network.dyadcount(y0)-network.edgecount(y0-is.na(y1))-network.naedgecount(y1)))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

y0<-network.initialize(n)
set.seed(321)
y0<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2))

cat("Complete data:\n")

set.seed(432)
y1<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2))

set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, estimate="CMLE", control=control.stergm(CMLE.control=control.ergm(MCMLE.conv.min.pval=0.8)))

stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

cat("Missing data:\n")

y1m<-network.copy(y1)
set.seed(654)
y1m[as.edgelist(simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))]<-NA

set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, estimate="CMLE", control=control.stergm(CMLE.control=control.ergm(MCMLE.conv.min.pval=0.8)))

stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))
