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

S<-5000

n<-20
m<-8
meanstats<-edges<-20
duration<-10
gamma<-logit(1-1/duration)

### Undirected

dyads<-n*(n-1)/2
density<-edges/dyads
theta<--log(((1+exp(gamma))/(density/(1-density)))-1)

cat("\nUndirected:\n")

g0<-network.initialize(n,dir=FALSE)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta)
print(gamma)

# Simulate from the fit.
dynsim<-simulatedyn(g1~edges,dissolve=~edges,dissolve.order="FormAndDiss",theta=theta,gamma=gamma,nsteps=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)

### Directed

dyads<-n*(n-1)
density<-edges/dyads
theta<--log(((1+exp(gamma))/(density/(1-density)))-1)

cat("\nDirected:\n")

g0<-network.initialize(n,dir=TRUE)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta)
print(gamma)

# Simulate from the fit.
dynsim<-simulatedyn(g1~edges,dissolve=~edges,dissolve.order="FormAndDiss",theta=theta,gamma=gamma,nsteps=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)

### Bipartite

dyads<-(n-m)*m
density<-edges/dyads
theta<--log(((1+exp(gamma))/(density/(1-density)))-1)

cat("\nBipartite:\n")

g0<-network.initialize(n,bipartite=m)

# Get a reasonably close starting network.
g1<-san(g0~edges,meanstats=meanstats,verbose=TRUE)

print(theta)
print(gamma)

# Simulate from the fit.
dynsim<-simulatedyn(g1~edges,dissolve=~edges,dissolve.order="FormAndDiss",theta=theta,gamma=gamma,nsteps=S,verbose=TRUE)

print.sim.stats(dynsim,meanstats,duration)


