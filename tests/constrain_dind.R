#  File tests/constrain_dind.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

library(ergm)

library(statnet.common)

mean_mat <- function(Mmin, Mmax){
  Mmin <- NVL(Mmin, Mmax)
  Mmax <- NVL(Mmax, Mmin)
  matrix(ifelse(rbinom(length(Mmin), 1, .5), Mmin, Mmax), nrow(Mmin), ncol(Mmin))
}

test_dind_constr <- function(y0, con, Mmin=NULL, Mmax=NULL, response=NULL, ...){
  nn <- network.dyadcount(y0, FALSE)
  if(!is.null(Mmin)){
    ymin <- simulate(NVL2(response, y0~sum, y0~edges), coef=-100, constraints=con, control=control.simulate.formula(MCMC.burnin=nn*100), response=response, ...)
    stopifnot(all(na.omit(c(suppressWarnings(as.matrix(ymin, attrname=response))==Mmin))))
  }
  if(!is.null(Mmax)){
    ymax <- simulate(NVL2(response, y0~sum, y0~edges), coef=+100, constraints=con, control=control.simulate.formula(MCMC.burnin=nn*100), response=response, ...)
    stopifnot(all(na.omit(c(as.matrix(ymax, attrname=response)==Mmax))))
  }
}

n <- 10
m <- 7

###### Unconstrained ######

Mmin <- matrix(0,n,n)
Mmax <- matrix(1,n,n)
diag(Mmax)<-0

#### Directed ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=TRUE)

test_dind_constr(y0, ~., Mmin, Mmax)

#### Undirected ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=FALSE)

test_dind_constr(y0, ~., Mmin, Mmax)

#### Unobserved ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=TRUE)
y0[2,3]<-NA
y0[2,10]<-NA

Mmin <- Mmax <- as.matrix(y0)

Mmin[2,10] <- Mmin[2,3] <- 0
Mmax[2,10] <- Mmax[2,3] <- 1

test_dind_constr(y0, ~observed, Mmin, Mmax)

#### Bipartite ####

Mmin <- matrix(0,m,n-m)
Mmax <- matrix(1,m,n-m)

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=FALSE, bipartite=m)

test_dind_constr(y0, ~., Mmin, Mmax)

#### Bipartite Unobserved ####

y0[7,8]<-NA
y0[6,9]<-NA

Mmin <- Mmax <- as.matrix(y0)

Mmin[6,2] <- Mmin[7,1] <- 0
Mmax[6,2] <- Mmax[7,1] <- 1

test_dind_constr(y0, ~observed, Mmin, Mmax) # in the block OR unobserved



######### Block diagonal #########

a <- rep(1:4,1:4)

Mmax <- Mmin <- matrix(0,n,n)

for(i in unique(a)){
  Mmax[a==i,a==i]<-1
}
diag(Mmax)<-0

#### Directed ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=TRUE)
y0 %v% "b" <- a

test_dind_constr(y0, ~blockdiag("b"), Mmin, Mmax)

#### Undirected ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=FALSE)
y0 %v% "b" <- a

test_dind_constr(y0, ~blockdiag("b"), Mmin, Mmax)

#### Unobserved: disjunction ####

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=TRUE)
y0 %v% "b" <- a
y0[2,3]<-NA
y0[2,10]<-NA

Mmin[2,10] <- Mmin[2,3] <- 0
Mmax[2,10] <- Mmax[2,3] <- 1

test_dind_constr(y0, ~blockdiag("b")-observed, Mmin, Mmax) # in the block OR unobserved

#### Unobserved: conjunction ####

Mmin <- Mmax <- as.matrix(y0)

Mmin[2,3] <- 0
Mmax[2,3] <- 1

test_dind_constr(y0, ~blockdiag("b")+observed, Mmin, Mmax) # in the block AND unobserved

#### Bipartite ####

a <- unlist(split(a, rep(1:2, n/2)))
a <- c(sort(a[1:m]), sort(a[-(1:m)]))

Mmin <- Mmax <- matrix(0,m,n-m)

for(i in unique(a)){
  Mmax[a[1:m]==i,a[-(1:m)]==i]<-1
}

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=FALSE, bipartite=m)
y0 %v% "b" <- a

test_dind_constr(y0, ~blockdiag("b"), Mmin, Mmax)

#### Bipartite Unobserved: disjunction ####

y0[7,8]<-NA
y0[6,9]<-NA

Mmin[6,2] <- Mmin[7,1] <- 0
Mmax[6,2] <- Mmax[7,1] <- 1

test_dind_constr(y0, ~blockdiag("b")-observed, Mmin, Mmax) # in the block OR unobserved

#### Bipartite Unobserved: conjunction ####

Mmin <- Mmax <- as.matrix(y0)

Mmin[6,2] <- 0
Mmax[6,2] <- 1

test_dind_constr(y0, ~blockdiag("b")+observed, Mmin, Mmax) # in the block AND unobserved


#### Valued ####

a <- rep(1:4,1:4)

Mmax <- Mmin <- matrix(0,n,n)

for(i in unique(a)){
  Mmax[a==i,a==i]<-4
}
diag(Mmax)<-0

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=TRUE, ignore.eval=FALSE, names.eval="w")
y0 %v% "b" <- a

test_dind_constr(y0, ~blockdiag("b"), Mmin, Mmax, response="w", reference=~DiscUnif(0,4))

