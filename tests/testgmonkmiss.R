library(ergm)
data(sampson)

#
# Create random 25% missing
#
# msamplike <- rergm(network.size(samplike),prob=0.25, directed=FALSE)
# samplike <- set.graph.attribute(samplike, "design", msamplike)
# summary(samplike)

respondent <- rep(FALSE,network.size(samplike))
respondent[sample(1:network.size(samplike), size=network.size(samplike)-2,replace=FALSE)] <- TRUE
respondent

summary(samplike)

degreedist(samplike)

efit <- ergm(samplike~edges + triangle, estimate="MPLE")
summary(efit)

efit <- ergm(samplike~edges + triangle, control=control.ergm(MCMLE.maxit=3,
      MCMC.samplesize=1000, MCMC.interval=1000, MCMC.burnin=1000))
summary(efit)

## Test bounded degrees.
efit <- ergm(samplike~edges + triangle, constraints=~bd(maxout=9))
summary(efit)

samplike <- set.vertex.attribute(samplike, "respondent", respondent)
rm(respondent)
summary(samplike)

efit <- ergm(samplike~edges + triangle, estimate="MPLE")
summary(efit)

efit <- ergm(samplike~edges + triangle, control=control.ergm(MCMLE.maxit=3,
    MCMC.samplesize=1000, MCMC.interval=1000, MCMC.burnin=1000))
summary(efit)

