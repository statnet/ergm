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

S<-100000

n<-200
m<-100
target.stats<-edges<-100
duration<-1000
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
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,nsim=S,verbose=TRUE)

#print.sim.stats(dynsim,target.stats,duration)

### Directed

dyads<-n*(n-1)
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)

cat("\nDirected:\n")

g0<-network.initialize(n,dir=TRUE)

# Get a reasonably close starting network.
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

print(coef.form)
print(coef.diss)

# Simulate from the fit.
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,nsim=S,verbose=TRUE)

#print.sim.stats(dynsim,target.stats,duration)

### Bipartite undirected

dyads<-(n-m)*m
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)

cat("\nBipartite:\n")

g0<-network.initialize(n,bipartite=m,directed=FALSE)

# Get a reasonably close starting network.
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

print(coef.form)
print(coef.diss)

# Simulate from the fit.
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,nsim=S,verbose=TRUE)

#print.sim.stats(dynsim,target.stats,duration)


