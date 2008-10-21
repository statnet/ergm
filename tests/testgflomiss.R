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

#efit <- ergm(flomarriage ~ edges + kstar(2), MCMCsamplesize=1000, interval=1000)
efit <- ergm(flomarriage ~ edges + kstar(2), MPLEonly=TRUE)
summary(efit)

flomarriage <- set.vertex.attribute(flomarriage, "respondent", respondent)
rm(respondent)
summary(flomarriage)

efit <- ergm(flomarriage ~ edges + kstar(2), MPLEonly=T)
summary(efit)

efit <- ergm(flomarriage ~ edges + kstar(2), MCMCsamplesize=1000, interval=1000)
summary(efit)

efit <- ergm(flomarriage ~ edges + kstar(2), theta=c(-1.6,0),startatMPLE=F)

#
# edges  -1.6     -1.74242 0.8557   0.044   0.041373 
# kstar2  0.0      0.02746 0.1774   0.877   0.009076 
#
summary(efit)
