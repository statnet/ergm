#  File tests/ergm.update.formula.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(ergm)

while(exists("test.network"))
  rm(test.network)

lhs.subst.san <- function(rhs,target.stats) {
        test.network <- network.initialize(n=10,directed=F)
        form <- ergm.update.formula(rhs,test.network~., from.new="test.network")
        san(form,target.stats=target.stats)
        }



stopifnot(summary(lhs.subst.san(~edges,20)~edges)==20)

# Test code created by Nicola Soriani
# depends on proper functioning of ergm.update.formula

predict.ergm.model<- function(model)  # only for Directed Network
{
 net<-model$network
 n<-network.size(net)
 beta<-model$coef
 eta <- matrix(0, n, n)


 for(i in 1:n)
	for(j in 1:n)
	{
	 alternative<-net
	 if(is.na(net[i,j]))
	  {
	   net[i,j]<-1
	   alternative[i,j]<-1 - net[i,j]
	   u2<- summary(ergm.update.formula(model$formula, alternative ~ ., from.new="alternative"))
	   delta<-u2-u1
	   prob<-1/(1+exp(sum(beta*delta)))
	   if (net[i,j]==1) eta[i,j]<- prob
       else eta[i,j]<-1-prob
       if (i==j) eta[i,j]<-0
       net[i,j]<-NA
       }
    else
     {
      alternative[i,j]<-1 - net[i,j]
      u1<-summary(ergm.update.formula(model$formula, net ~ ., from.new="net"))
      u2<-summary(ergm.update.formula(model$formula, alternative ~ ., from.new="alternative"))
      delta<- u2-u1
      prob<- 1/(1+exp(sum(beta*delta)))
      if (net[i,j]==1) eta[i,j]<- prob
      else eta[i,j]<-1-prob
      if (i==j) eta[i,j]<-0
      }
    }
  return(eta)
}


opttest({
#data(florentine)
data(g4)

#Ergm<-ergm(flobusiness~edges+gwesp(1,fixed=TRUE))
Ergm<-ergm(g4~edges+gwesp(1,fixed=TRUE))

stopifnot(all(diag(predict.ergm.model(Ergm))==0))

}, "part of ergm.update.formula")
