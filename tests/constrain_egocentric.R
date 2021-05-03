#  File tests/constrain_blockdiag.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

library(ergm)
n <- 10
m <- 7

propw <- "TNT"
#propw <- "random"

repeat{
  a <- as.logical(rbinom(n, 1, .5))
  Mb <- Mi <- Mo <- matrix(0L,n,n)
  
  Mo[!a,] <- 1L
  Mi[,!a] <- 1L
  Mb <- Mo * Mi
  
  diag(Mo) <- diag(Mi) <- diag(Mb) <-0L

  if(any(Mo) && any(Mi) && any(Mb)) break
}
#### Directed ####

y0 <- network.initialize(n, directed=TRUE)
y0 %v% "a" <- a

y <- simulate(y0~edges, coef=100, constraints=~egocentric("a", dir="o"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))
stopifnot(all(as.matrix(y)==Mo))

y <- simulate(y0~edges, coef=100, constraints=~egocentric("a", dir="i"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))
stopifnot(all(as.matrix(y)==Mi))

y <- simulate(y0~edges, coef=100, constraints=~egocentric("a"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))
stopifnot(all(as.matrix(y)==Mb))


#### Undirected ####

y0 <- network.initialize(n, directed=FALSE)
y0 %v% "a" <- a

y <- simulate(y0~edges, coef=100, constraints=~egocentric("a"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))
stopifnot(all(as.matrix(y)==Mb))

## #### Unobserved ####

## y0 <- network.initialize(n, directed=TRUE)
## y0 %v% "a" <- a
## y0[2,3]<-NA
## y0[2,10]<-NA

## y <- simulate(y0~edges, coef=100, constraints=~egocentric("a")+observed, control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))

## M[]<-0
## M[2,3]<-1

## stopifnot(all(as.matrix(y)==M))

#### Bipartite ####

y0 <- network.initialize(n, directed=FALSE, bipartite=m)

repeat{
  a <- as.logical(rbinom(n, 1L, .5))
  ae <- a[seq_len(m)]
  aa <- a[m+seq_len(n-m)]
  
  M <- matrix(1L,m,n-m)
  
  M[ae,] <- 0L
  M[,aa] <- 0L
  
  if(any(M)) break
}

y0 %v% "a" <- a

y <- simulate(y0~edges, coef=100, constraints=~egocentric("a"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))
stopifnot(all(as.matrix(y)==M))

## #### Bipartite Unobserved ####

## y0 <- network.initialize(n, directed=FALSE, bipartite=m)
## y0 %v% "b" <- a
## y0[7,8]<-NA
## y0[6,9]<-NA

## y <- simulate(y0~edges, coef=100, constraints=~blockdiag("b")+observed, control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))

## M[]<-0
## M[6,2]<-1

## stopifnot(all(as.matrix(y)==M))

## #### Multiple ####

## n <- 10
## a1 <- rep(1:4,1:4)
## a2 <- rep(1:2,each=5)

## M1<- matrix(0,n,n)
## for(i in unique(a1)){
##   M1[a1==i,a1==i]<-1
## }
## diag(M1)<-0

## M2<- matrix(0,n,n)
## for(i in unique(a2)){
##   M2[a2==i,a2==i]<-1
## }
## diag(M2)<-0

## M <- M1*M2

## y0 <- network.initialize(n, directed=FALSE)
## y0 %v% "b1" <- a1
## y0 %v% "b2" <- a2

## y <- simulate(y0~edges, coef=100, constraints=~blockdiag("b1") + blockdiag("b2"), control=control.simulate.formula(MCMC.burnin=10000,MCMC.prop.weight=propw))

## stopifnot(all(as.matrix(y)==M))

