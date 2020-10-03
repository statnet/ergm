#  File tests/constrain_dind.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
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

#### Dyads operator ####
## TODO: Put in the same framework as the others.

data(sampson)
fix_g <- coef(ergm(samplike~edges, constraints=~Dyads(~nodematch("group"))))
vary_g <- coef(ergm(samplike~edges, constraints=~Dyads(vary=~nodematch("group"))))
fix_g_and_c <- coef(ergm(samplike~edges, constraints=~Dyads(~nodematch("group")+nodematch("cloisterville"))))
fix_g_vary_c <- coef(ergm(samplike~edges, constraints=~Dyads(~nodematch("group"),~nodematch("cloisterville"))))
vary_g_or_c <- coef(ergm(samplike~edges, constraints=~Dyads(vary=~nodematch("group")+nodematch("cloisterville"))))
vary_g_fix_c <- coef(ergm(samplike~edges, constraints=~Dyads(vary=~nodematch("group"))+Dyads(~nodematch("cloisterville"))))

# Check:
m <- as.matrix(samplike)
g <- outer(samplike%v%"group",samplike%v%"group",FUN=`==`)
c <- outer(samplike%v%"cloisterville",samplike%v%"cloisterville",FUN=`==`)
n <- network.size(samplike)
logit <- function(p) log(p/(1-p))
stopifnot(isTRUE(all.equal(fix_g,logit(sum((!g)*m)/(sum(!g))),check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(vary_g,logit(sum(g*m)/(sum(g)-n)),check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(fix_g_and_c,logit(sum((!g&!c)*m)/(sum(!g&!c))),check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(fix_g_vary_c,logit(sum((!g|c)*m)/(sum(!g|c)-n)),check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(vary_g_or_c,logit(sum((g|c)*m)/(sum(g|c)-n)),check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(vary_g_fix_c,logit(sum((g&!c)*m)/(sum(g&!c))),check.attributes=FALSE)))
