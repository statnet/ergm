# Test code created by Nicola Soriani
# depends on proper functioning of update.formula

library(ergm)
predict.ergm.model<- function(model)  # only for Directed Network
{
 net<-model$"network"
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
	   u2<- summary(update.formula(model$formula, alternative ~ .))
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
      u1<-summary(update.formula(model$formula, net ~ .))
      u2<-summary(update.formula(model$formula, alternative ~ .))
      delta<- u2-u1
      prob<- 1/(1+exp(sum(beta*delta)))
      if (net[i,j]==1) eta[i,j]<- prob
      else eta[i,j]<-1-prob
      if (i==j) eta[i,j]<-0
      }
    }
  return(eta)
}

data(florentine)

Ergm<-ergm(flobusiness~edges+gwesp(1,fixed=TRUE))

if (!all(diag(predict.ergm.model(Ergm))==0)) stop("failed update.formula test")

