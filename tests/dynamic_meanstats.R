library(ergm)
library(coda)
#ergm<-ergm2

n<-50

g0<-network(n,dir=FALSE)

#            meandeg, degree(1)
meanstats<-c(      1,    n*0.6)

# Get a reasonably close starting network.
g1<-san(g0~meandeg+degree(1),meanstats=meanstats,verbose=TRUE)

# Fit the model.
dynfit<-ergm2(g1~meandeg+degree(1),dissolve=g1~edges,gamma=log(.95/.05),meanstats=meanstats,control=ergm.control(style="Robbins-Monro"),verbose=TRUE)

theta<-dynfit$coef
print(theta)

gamma<-log(.95/.05)
print(gamma)

# Simulate from the fit.
dynsim<-simulatedyn(g1~meandeg+degree(1),dissolve=g1~edges,theta=theta,gamma=gamma,nsteps=1000,verbose=TRUE)

dynsim.gf<-ergm.godfather(g1~meandeg+degree(1),sim=dynsim,verbose=TRUE)

# Compare each time point's network statistics as returned by
# simulation and as returned by a "replay" of the simulation using the
# Godfather Proposal.
# If they don't match, it's a bug.
print(all.equal(dynsim$stats.form,dynsim.gf$stats))

# Print out the resulting meanstats and their t-value w.r.t. the
# target meanstats.
print(meanstats)
print(apply(dynsim$stats.form,2,mean))
print((apply(dynsim$stats.form,2,mean)-meanstats)/sqrt(apply(dynsim$stats.form,2,var)/effectiveSize(mcmc(dynsim$stats.form))))
print(apply(dynsim.gf$stats,2,mean))
print((apply(dynsim.gf$stats,2,mean)-meanstats)/sqrt(apply(dynsim.gf$stats,2,var)/effectiveSize(mcmc(dynsim.gf$stats))))

print(mean(duration.matrix(dynsim)$duration))


