library(ergm)
n<-20
base.net <- network.initialize(n=n,directed=FALSE)
n*c(1,.2,.4)
ergm.fit<-ergm(base.net~edges+degree(c(0,1)),meanstats=n*c(1,.2,.4))
summary(ergm.fit)
