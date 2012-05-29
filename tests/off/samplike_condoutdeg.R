library(ergm)
data(sampson)


degreedist(samplike)

outdegrees <- apply(as.matrix(samplike, m="a"), 1, sum)
table(outdegrees)

efit <- ergm(samplike ~ edges + triangle, estimate="MPLE")
summary(efit)

#
# This fit holds the out degrees fixed
#
efit <- ergm(samplike ~ edges + triangle, constraints=~outdegrees,
  control=control.ergm(MCMLE.maxit=3, MCMC.samplesize=10000))
summary(efit)

