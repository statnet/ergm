library(ergm)

n <- 50

logit<-function(p) log(p/(1-p))

y <- network.initialize(n, directed=FALSE) # Create an empty network
y <- simulate(y~edges, coef=logit(0.12), control=control.simulate(MCMC.burnin=2*n^2))
y.miss <- simulate(y~edges, coef=logit(0.01))
y[as.edgelist(y.miss)] <- NA

cdfit<-ergm(y~edges+gwesp(), estimate="CD", control=control.ergm(CD.nsteps=50, MCMC.samplesize=100))
