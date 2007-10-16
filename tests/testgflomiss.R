library(ergm)
#
data(gflo)

##
## Create random 5% missing
##
#mgflo <- rergm(graph.size(gflo),prob=0.05, directed=FALSE)
#gflo <- set.graph.attribute(gflo, "design", mgflo)
#summary(gflo)

#
# Create random 2 nodes who are non-respondents
#
#respondent <- rmultinom(n=1, size=graph.size(gflo)-2,
#                        prob=rep(1,graph.size(gflo)))
#respondent

respondent <- rep(FALSE,graph.size(gflo))
respondent[sample(1:graph.size(gflo), size=graph.size(gflo)-2,replace=FALSE)] <- TRUE
respondent

#
#one <- matrix(1,ncol=1,nrow=graph.size(gflo))
#orespondent <- one %*% t(respondent) + respondent %*% t(one) - respondent %*% t(respondent)
#orespondent <- 1-orespondent
#diag(orespondent) <- 0
##
#mgflo <- graph(orespondent, directed=FALSE)
#summary(mgflo)
#sociomatrix(mgflo)
#gflo <- set.graph.attribute(gflo, "design", mgflo)

#efit <- ergm(gflo ~ edges + kstar(2), MCMCsamplesize=1000, interval=1000)
efit <- ergm(gflo ~ edges + kstar(2), MPLEonly=TRUE)
summary(efit)

gflo <- set.vertex.attribute(gflo, "respondent", respondent)
rm(respondent)
summary(gflo)

efit <- ergm(gflo ~ edges + kstar(2), MPLEonly=T)
summary(efit)

efit <- ergm(gflo ~ edges + kstar(2), MCMCsamplesize=1000, interval=1000)
summary(efit)

efit <- ergm(gflo ~ edges + kstar(2), theta=c(-1.6,0),startatMPLE=F)

#
# edges  -1.6     -1.74242 0.8557   0.044   0.041373 
# kstar2  0.0      0.02746 0.1774   0.877   0.009076 
#
summary(efit)
