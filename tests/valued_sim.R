library(ergm)

## Poisson-reference
cat("Poisson-reference ERGM\n")
load("testnet3u.RData")

theta<-1
cat("Target mean:",exp(theta),"\n",sep="")

s<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE)

cat("Simulated mean (statsonly):",mean(s)/3,"\n",sep="")

s.full<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE)

cat("Simulated mean (full, computed):",mean(sapply(s.full,function(x) sum(x%e%"w")))/3,"\n",sep="")
cat("Simulated mean (full, stats):",mean(attr(s.full,"stats"))/3,"\n",sep="")

## Poisson-reference, zero-inflated
cat("Poisson-reference ERGM with zero-inflation\n")

theta<-1
cat("Target mean:",exp(theta),"\n",sep="")

s<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE,
            control=control.simulate(prop.weights="0inflated"))

cat("Simulated mean (statsonly):",mean(s)/3,"\n",sep="")

s.full<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE,
                 control=control.simulate(prop.weights="0inflated"))

cat("Simulated mean (full, computed):",mean(sapply(s.full,function(x) sum(x%e%"w")))/3,"\n",sep="")
cat("Simulated mean (full, stats):",mean(attr(s.full,"stats"))/3,"\n",sep="")


## StdNormal-reference
cat("Standard-normal-reference ERGM with mutuality by correlation\n")
load("testnet3d.RData")

# Joint distribution of (i,j) and (j,i)
mu<-1
sig<-2
rho<-.3

# Naural parameters of bivariate normal
# We ought to be able to specify them in a more interpretable way.
denom<--2*(1-rho^2)*sig^2
xx.coef<-1/denom+1/2 # 1/2 is from the reference measure
x.coef<-2*(-1+rho)*mu/denom
xy.coef<--2*rho/denom

theta<-c(x.coef,xy.coef,xx.coef)

cat("mean=",mu,", var=",sig^2,", corr=",rho,"\neta=(",paste(theta,collapse=","),")\n",sep="")

s<-simulate(testnet3d~sum+mutual("product")+sum(pow=2),reference="StdNormal",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE)

cat("Simulated mean (statsonly):",mean(s[,1])/6,"\n",sep="")

s.full<-simulate(testnet3d~sum+mutual("product")+sum(pow=2),reference="StdNormal",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE)

s.cells<-sapply(s.full, function(x) as.matrix(x,m="a",a="w"),simplify=FALSE)
cat("Simulated means (target=",mu,"):\n",sep="")
print(matrix(c(NA,
  mean(sapply(s.cells,"[",1,2)),
  mean(sapply(s.cells,"[",1,3)),
  mean(sapply(s.cells,"[",2,1)),
  NA,
  mean(sapply(s.cells,"[",2,3)),
  mean(sapply(s.cells,"[",3,1)),
  mean(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated vars (target=",sig^2,"):\n",sep="")
print(matrix(c(NA,
  var(sapply(s.cells,"[",1,2)),
  var(sapply(s.cells,"[",1,3)),
  var(sapply(s.cells,"[",2,1)),
  NA,
  var(sapply(s.cells,"[",2,3)),
  var(sapply(s.cells,"[",3,1)),
  var(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated correlations (1,2) (1,3) (2,3) (target=",rho,"):\n",sep="")
print(c(cor(sapply(s.cells,"[",1,2),sapply(s.cells,"[",2,1)),
        cor(sapply(s.cells,"[",1,3),sapply(s.cells,"[",3,1)),
        cor(sapply(s.cells,"[",2,3),sapply(s.cells,"[",3,2))))

## StdNormal-reference with rank constraint
cat("Standard-normal-reference ERGM with rank constraint\n")
load("testrank3d.RData")

s.full<-simulate(testrank3d~sum,reference="StdNormal",response="w",theta0=0,burnin=10000,nsim=1000,statsonly=FALSE,constraints=~ranks)
s.cells<-sapply(s.full, function(x) as.matrix(x,m="a",a="w"),simplify=FALSE)
ref.sample<-pmax(rnorm(10000),rnorm(10000))

cat("Simulated means (target[1:2,]=+-",mean(ref.sample),";target[3,]=0):\n",sep="")
print(matrix(c(NA,
  mean(sapply(s.cells,"[",1,2)),
  mean(sapply(s.cells,"[",1,3)),
  mean(sapply(s.cells,"[",2,1)),
  NA,
  mean(sapply(s.cells,"[",2,3)),
  mean(sapply(s.cells,"[",3,1)),
  mean(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))

cat("Simulated vars (target[1:2,]=+-",var(ref.sample),";target[3,]=0):\n",sep="")
print(matrix(c(NA,
  var(sapply(s.cells,"[",1,2)),
  var(sapply(s.cells,"[",1,3)),
  var(sapply(s.cells,"[",2,1)),
  NA,
  var(sapply(s.cells,"[",2,3)),
  var(sapply(s.cells,"[",3,1)),
  var(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))
