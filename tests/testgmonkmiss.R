library(ergm)
data(gmonk)

#
# Create random 25% missing
#
# mgmonk <- rergm(graph.size(gmonk),prob=0.25, directed=FALSE)
# gmonk <- set.graph.attribute(gmonk, "design", mgmonk)
# summary(gmonk)

respondent <- rep(FALSE,graph.size(gmonk))
respondent[sample(1:graph.size(gmonk), size=graph.size(gmonk)-2,replace=FALSE)] <- TRUE
respondent

summary(gmonk)

degreedist(gmonk)

efit <- ergm(gmonk~edges + triangle, MPLEonly=T)
summary(efit)

efit <- ergm(gmonk~edges + triangle, maxit=3,
 MCMCsamplesize=1000, interval=1000, burnin=1000)
summary(efit)

gmonk <- set.vertex.attribute(gmonk, "respondent", respondent)
rm(respondent)
summary(gmonk)

efit <- ergm(gmonk~edges + triangle, MPLEonly=T)
summary(efit)

efit <- ergm(gmonk~edges + triangle, maxit=3,
 MCMCsamplesize=1000, interval=1000, burnin=1000)
summary(efit)
