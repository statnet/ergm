load("testnet3u.RData")
load("testrank3d.RData")
library(ergm)

## Poisson-reference

theta<-1

print(exp(theta))

s<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE)

print(mean(s)/3) # Should be around exp(theta).

s.full<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE)

print(mean(sapply(s.full$networks,function(x) sum(x%e%"w")))/3) # Should also be around exp(theta)
print(mean(s.full$stats)/3)  # Should be equal to the previous line.

## StdNormal-reference

theta<--1

print(theta)

s<-simulate(testnet3u~sum,reference="StdNormal",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE)

print(mean(s)/3) # Should be around theta.
print(sd(s)/sqrt(3)) # Should be around 1.

s.full<-simulate(testnet3u~sum,reference="StdNormal",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE)

print(mean(sapply(s.full$networks,function(x) sum(x%e%"w")))/3) # Should also be around theta
print(mean(s.full$stats)/3)  # Should be equal to the previous line.
print(sd(sapply(s.full$networks,function(x) sum(x%e%"w")))/sqrt(3)) # Should also be around 1
print(sd(s.full$stats)/sqrt(3)) # Should be equal to the previous line.

## StdNormal-reference with rank constraint

s.full<-simulate(testrank3d~sum,reference="StdNormal",response="w",theta0=0,burnin=10000,nsim=1000,statsonly=FALSE,constraints=~ranks)
s.cells<-sapply(s.full$networks, function(x) as.matrix(x,m="a",a="w"),simplify=FALSE)
ref.sample<-pmax(rnorm(10000),rnorm(10000))
print(mean(ref.sample))
print(sd(ref.sample))
print(matrix(c(NA,
  mean(sapply(s.cells,"[",1,2)),
  mean(sapply(s.cells,"[",1,3)),
  mean(sapply(s.cells,"[",2,1)),
  NA,
  mean(sapply(s.cells,"[",2,3)),
  mean(sapply(s.cells,"[",3,1)),
  mean(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))
print(matrix(c(NA,
  sd(sapply(s.cells,"[",1,2)),
  sd(sapply(s.cells,"[",1,3)),
  sd(sapply(s.cells,"[",2,1)),
  NA,
  sd(sapply(s.cells,"[",2,3)),
  sd(sapply(s.cells,"[",3,1)),
  sd(sapply(s.cells,"[",3,2)),
  NA),3,3,byrow=TRUE))
