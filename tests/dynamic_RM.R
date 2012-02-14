library(ergm)
library(coda)

n<-50

g0<-network.initialize(n,dir=FALSE)

#            meandeg, degree(1)
target.stats<-c(      n*1/2,    n*0.6)

# Get a reasonably close starting network.
g1<-san(g0~meandeg+degree(1),target.stats=target.stats,verbose=TRUE)

# Fit the model.
dynfit<-stergm(g1,formation=~edges+degree(1),dissolution=~offset(edges), targets="formation", estimate="EGMoME", offset.coef.diss=log(.95/.05),target.stats=target.stats,verbose=TRUE,control=control.stergm(RM.interval=100))

coef.form<-dynfit$formation.fit$coef
print(coef.form)

coef.diss<-log(.95/.05)
print(coef.diss)

# Simulate from the fit.
dynsim<-simulate(dynfit,nsim=1000,verbose=TRUE)

dynsim.gf<-ergm.godfather(g1~edges+degree(1),sim=dynsim,verbose=TRUE)

# Compare each time point's network statistics as returned by
# simulation and as returned by a "replay" of the simulation using the
# Godfather Proposal.
# If they don't match, it's a bug.
# Note that the first row of the Godfather stats is the initial network.
if(!isTRUE(all.equal(dynsim$stats.form,dynsim.gf$stats[-1,])))
  stop("Formation statistics returned by simulate differ from those returned",
       "by ergm.godfather. This is a bug.")

# Print out the resulting target.stats and their t-value w.r.t. the
# target target.stats.
print(target.stats)
target.stats.sim<-apply(dynsim$stats.form,2,mean)
print(target.stats.sim)
print(effectiveSize(mcmc(dynsim$stats.form)))
print((target.stats.sim-target.stats)/sqrt(apply(dynsim$stats.form,2,var)/effectiveSize(mcmc(dynsim$stats.form))))

print(mean(duration.matrix(dynsim)$duration))

# Simulate from an equivalent fit.
dynsim<-simulate(g1,~edges+degree(1),dissolution=~dyadcov(matrix(1,n,n))+edges,coef.form=coef.form,coef.diss=c(1,coef.diss-1),nsim=1000,verbose=TRUE)

print(mean(duration.matrix(dynsim)$duration))

