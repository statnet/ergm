#  File ergm/tests/dynamic_simtest_undirected.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
library(coda)

logit<-function(p)log(p/(1-p))

# NB:  duration.matrix function no longer exists in ergm package
#print.sim.stats<-function(dynsim,m,d){
#  t.score<-function(x,m) (mean(x)-m)/sqrt(apply(cbind(x),2,var)/effectiveSize(mcmc(x)))
#  target.stats.sim<-apply(dynsim$stats.form,2,mean)
#  durations<-duration.matrix(dynsim)$duration
#  cat('Edge count:\n   Target:',m,', Simulated:',target.stats.sim,', t:', t.score(dynsim$stats.form,m) ,'\n')
#  cat('Duration:\n   Target:',d,', Simulated:',mean(durations),', t:', t.score(durations,d) ,'\n')
#}

coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)

S<-1000

n<-200
m<-100
target.stats<-edges<-100
duration<-100
coef.diss<-logit(1-1/duration)

### Undirected

dyads<-n*(n-1)/2
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)

cat("\nUndirected:\n")

g0<-network.initialize(n,dir=FALSE)

# Get a reasonably close starting network.
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

print(coef.form)
print(coef.diss)

# Simulate from the fit.
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE)

#print.sim.stats(dynsim,target.stats,duration)
