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
dynfit<-ergm2(g1~meandeg+degree(1),dissolve=g1~edges,dissolve.orger="FormAndDiss",gamma=log(.95/.05),meanstats=meanstats,control=ergm.control(style="Robbins-Monro"),verbose=2)

theta<-dynfit$coef
print(theta)

gamma<-log(.95/.05)
print(gamma)

# Simulate from the fit.
dynsim<-simulatedyn(g1~meandeg+degree(1),dissolve=g1~edges,dissolve.order="FormAndDiss",theta=theta,gamma=gamma,nsteps=1000,verbose=TRUE)

dynsim.gf<-ergm.godfather(g1~meandeg+degree(1),sim=dynsim,verbose=TRUE)

# Compare each time point's network statistics as returned by
# simulation and as returned by a "replay" of the simulation using the
# Godfather Proposal.
# If they don't match, it's a bug.
# Note that the first row of the Godfather stats is the initial network.
if(!isTRUE(all.equal(dynsim$stats.form,dynsim.gf$stats[-1,])))
  stop("Formation statistics returned by simulatedyn differ from those returned",
       "by ergm.godfather. This is a bug.")

# Print out the resulting meanstats and their t-value w.r.t. the
# target meanstats.
print(meanstats)
meanstats.sim<-apply(dynsim$stats.form,2,mean)
print(meanstats.sim)
print((meanstats.sim-meanstats)/sqrt(apply(dynsim$stats.form,2,var)/effectiveSize(mcmc(dynsim$stats.form))))

print(mean(duration.matrix(dynsim)$duration))


