library(ergm) 
#
data(g4)

g4 <- 1*((g4 + t(g4)) > 0)

g4

mg4 <- g4 - g4
mg4[2,3] <- 1
mg4[3,2] <- 1
#mg4[c(1,2)] <- 1
#mg4[c(2,1)] <- 1

g4 <- graph(g4, directed=FALSE)

mg4 <- graph(mg4, directed=FALSE)
sociomatrix(mg4)
g4 <- set.graph.attribute(g4, "design", mg4)

g4

degreedist(g4)

#gout <- ergm(g4 ~ triangle ,
#             MPLEonly=TRUE)
#summary(gout)

#load(file="goutg4.RData")

gout <- ergm(g4 ~ triangle, MCMCsamplesize=1000, interval=1000, maxit=3)
summary(gout)

#save(gout, file="goutg4.RData")

rout <- rergm(gout)
summary(rout)
rout

rout <- rergm(g4 ~ triangle, MCMCsamplesize=1000, interval=1000,
 theta0=gout$coef)
summary(rout)
rout

