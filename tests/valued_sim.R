load("testnet3u.RData")
library(ergm)

## Poisson-reference

theta<-1

print(exp(theta))

s<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=TRUE)

print(mean(s)/3) # Should be around exp(theta).

s.full<-simulate(testnet3u~sum,reference="Poisson",response="w",theta0=theta,burnin=10000,nsim=1000,statsonly=FALSE)

print(mean(sapply(s.full$networks,function(x) sum(x%e%"w")))/3) # Should also be around exp(theta)
print(mean(s.full$stats)/3)  # Should be equal to the previous line.
