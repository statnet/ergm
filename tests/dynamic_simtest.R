library(ergm)
library(coda)

logit<-function(p)log(p/(1-p))

print.sim.stats<-function(dynsim,m,d){
  t.score<-function(x,m) (mean(x)-m)/sqrt(apply(cbind(x),2,var)/effectiveSize(mcmc(x)))
  meanstats.sim<-apply(dynsim$stats.form,2,mean)
  durations<-duration.matrix(dynsim)$duration
  cat('Edge count:\n   Target:',m,', Simulated:',meanstats.sim,', t:', t.score(dynsim$stats.form,m) ,'\n')
  cat('Duration:\n   Target:',d,', Simulated:',mean(durations),', t:', t.score(durations,d) ,'\n')
}

theta.form.f<-function(theta.diss,density) -log(((1+exp(theta.diss))/(density/(1-density)))-1)

S<-100000

n<-200
m<-100
meanstats<-edges<-100
duration<-1000
theta.diss<-logit(1-1/duration)

### Undirected

dyads<-n*(n-1)/2
density<-edges/dyads
theta.form<-theta.form.f(theta.diss,density)

cat("\nUndirected:\n")

g0<-network.initialize(n,dir=FALSE)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta.form)
print(theta.diss)

# Simulate from the fit.
dynsim<-simulate(g1~edges,dissolution=~edges,stergm.order="FormAndDiss",theta.form=theta.form,theta.diss=theta.diss,nsim=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)

### Directed

dyads<-n*(n-1)
density<-edges/dyads
theta.form<-theta.form.f(theta.diss,density)

cat("\nDirected:\n")

g0<-network.initialize(n,dir=TRUE)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta.form)
print(theta.diss)

# Simulate from the fit.
dynsim<-simulate(g1~edges,dissolution=~edges,stergm.order="FormAndDiss",theta.form=theta.form,theta.diss=theta.diss,nsim=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)

### Bipartite undirected

dyads<-(n-m)*m
density<-edges/dyads
theta.form<-theta.form.f(theta.diss,density)

cat("\nBipartite:\n")

g0<-network.initialize(n,bipartite=m,directed=FALSE)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta.form)
print(theta.diss)

# Simulate from the fit.
dynsim<-simulate(g1~edges,dissolution=~edges,stergm.order="FormAndDiss",theta.form=theta.form,theta.diss=theta.diss,nsim=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)


