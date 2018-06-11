#  File tests/testgflomiss.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
#
data(florentine)

##
## Create random 5% missing
##
#mflomarriage <- rergm(network.size(flomarriage),prob=0.05, directed=FALSE)
#flomarriage <- set.graph.attribute(flomarriage, "design", mflomarriage)
#summary(flomarriage)

#
# Create random 2 nodes who are non-respondents
#
#respondent <- rmultinom(n=1, size=network.size(flomarriage)-2,
#                        prob=rep(1,network.size(flomarriage)))
#respondent

respondent <- rep(FALSE,network.size(flomarriage))
respondent[sample(1:network.size(flomarriage), size=network.size(flomarriage)-2,replace=FALSE)] <- TRUE
respondent

#
#one <- matrix(1,ncol=1,nrow=network.size(flomarriage))
#orespondent <- one %*% t(respondent) + respondent %*% t(one) - respondent %*% t(respondent)
#orespondent <- 1-orespondent
#diag(orespondent) <- 0
##
#mflomarriage <- graph(orespondent, directed=FALSE)
#summary(mflomarriage)
#sociomatrix(mflomarriage)
#flomarriage <- set.graph.attribute(flomarriage, "design", mflomarriage)

#efit <- ergm(flomarriage ~ edges + kstar(2), control=control.ergm(MCMC.samplesize=1000, MCMC.interval=1000))
efit <- ergm(flomarriage ~ edges + kstar(2), estimate="MPLE")
summary(efit)

flomarriage <- set.vertex.attribute(flomarriage, "respondent", respondent)
rm(respondent)
summary(flomarriage)

efit <- ergm(flomarriage ~ edges + kstar(2), estimate="MPLE")
summary(efit)

efit <- ergm(flomarriage ~ edges + kstar(2), control=control.ergm(MCMC.samplesize=1000, MCMC.interval=1000))
summary(efit)

efit <- ergm(flomarriage ~ edges + kstar(2), control=control.ergm(init=c(-1.6,0)))

#
# edges  -1.6     -1.74242 0.8557   0.044   0.041373 
# kstar2  0.0      0.02746 0.1774   0.877   0.009076 
#
summary(efit)
},,"undirected network with missing data and dyadic dependence")
