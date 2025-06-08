#  File tests/testthat/test-constrain-blockdiag.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

unloadNamespace("ergm.count")

mean_mat <- function(Mmin, Mmax){
  Mmin <- statnet.common::NVL(Mmin, Mmax)
  Mmax <- statnet.common::NVL(Mmax, Mmin)
  matrix(ifelse(rbinom(length(Mmin), 1, .5), Mmin, Mmax), nrow(Mmin), ncol(Mmin))
}

test_dind_constr <- function(y0, con, Mmin=NULL, Mmax=NULL, response=NULL, ...){
  nn <- network.dyadcount(y0, FALSE)
  test_that(paste0("blockdiag constraint with constraint = ", format(con), ", and ", if(is.directed(y0)) "directed " else "undirected ", if(is.bipartite(y0)) "bipartite ", if(!is.null(response)) "valued ", "network"), {
    ymin <- simulate(statnet.common::NVL2(response, y0~sum, y0~edges), coef=-100, constraints=con, control=control.simulate.formula(MCMC.burnin=nn*100), response=response, ...)
    expect_true(all(na.omit(c(suppressWarnings(as.matrix(ymin, attrname=response))==Mmin))))

    ymax <- simulate(statnet.common::NVL2(response, y0~sum, y0~edges), coef=+100, constraints=con, control=control.simulate.formula(MCMC.burnin=nn*100), response=response, ...)
    expect_true(all(na.omit(c(as.matrix(ymax, attrname=response)==Mmax))))
  })
}

n <- 10
m <- 7

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

#### Unobserved: disjunction ergmlhs specification ####

y0l <- y0
y0l %ergmlhs% "constraints" <- ~blockdiag("b")
test_dind_constr(y0l, ~-observed, Mmin, Mmax) # dot by default
test_dind_constr(y0l, ~.-observed, Mmin, Mmax) # dot specified

#### Unobserved: conjunction ####

Mmin <- Mmax <- as.matrix(y0)

Mmin[2,3] <- 0
Mmax[2,3] <- 1

test_dind_constr(y0, ~blockdiag("b")+observed, Mmin, Mmax) # in the block AND unobserved

#### Unobserved: conjunction ergmlhs specification ####

test_dind_constr(y0l, ~+observed, Mmin, Mmax) # dot by default
test_dind_constr(y0l, ~.+observed, Mmin, Mmax) # dot specified

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

#### Multiple ####

a1 <- rep(1:4,1:4)
a2 <- rep(1:2,each=5)

Mmin <- M1 <- M2 <- matrix(0,n,n)
for(i in unique(a1)){
  M1[a1==i,a1==i]<-1
}
diag(M1)<-0

for(i in unique(a2)){
  M2[a2==i,a2==i]<-1
}
diag(M2)<-0

Mmax <- M1*M2

y0 <- as.network(mean_mat(Mmin,Mmax), matrix.type="adjacency", directed=FALSE)
y0 %v% "b1" <- a1
y0 %v% "b2" <- a2

test_dind_constr(y0, ~blockdiag("b1") + blockdiag("b2"), Mmin, Mmax)

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

library(ergm.count)
