context("test-valued-terms.R")

# Correct values. Note that for undirected networks, this needs to be
# divied by 2.
transitiveweights <- function(m,ties.f=pmin,combine.f=max,compare.f=min) sum(unlist(sapply(1:nrow(m),function(i) sapply(1:nrow(m),function(j) if(j==i) 0 else compare.f(m[i,j],combine.f(ties.f(m[i,-c(j,i)], m[-c(i,j),j])))))))
cyclicalweights <- function(m,ties.f=pmin,combine.f=max,compare.f=min) sum(unlist(sapply(1:nrow(m),function(i) sapply(1:nrow(m),function(j) if(j==i) 0 else compare.f(m[j,i],combine.f(ties.f(m[i,-c(j,i)], m[-c(i,j),j])))))))

pgeomean <- function(x,y) sqrt(x*y)

summary.call <- function(y)
  summary(y
          ~ transitiveweights("min","max","min")
          + transitiveweights("min","sum","min")
          + transitiveweights("geomean","sum","geomean")
          + cyclicalweights("min","max","min")
          + cyclicalweights("min","sum","min")
          + cyclicalweights("geomean","sum","geomean"),
          response="w")

simulate.call <- function(y)
  simulate(y
           ~ transitiveweights("min","max","min")
           + transitiveweights("min","sum","min")
           + transitiveweights("geomean","sum","geomean")
           + cyclicalweights("min","max","min")
           + cyclicalweights("min","sum","min")
           + cyclicalweights("geomean","sum","geomean"),
           coef=rep(0,6),reference=~DiscUnif(0,4),response="w",control=control.simulate(MCMC.burnin=0,MCMC.interval=100),nsim=100,output="pending_update_network")

test_that("valued triadic effects in undirected networks", {

  # Undirected
  y <- network.initialize(20, dir=FALSE)
  
  # Check the s_ statistics
  
  y <- simulate(y~sum,coef=0,reference=~DiscUnif(0,4),response="w",control=control.simulate(MCMC.burnin=1000),nsim=1)
  
  y.summ <- summary.call(y)
  
  y.m <- as.matrix(y, a="w")

  expect_equivalent(y.summ[1], transitiveweights(y.m, pmin, max, min)/2)
  expect_equivalent(y.summ[2], transitiveweights(y.m, pmin, sum, min)/2)
  expect_equivalent(y.summ[3], transitiveweights(y.m, pgeomean, sum, pgeomean)/2)
  expect_equivalent(y.summ[4], cyclicalweights(y.m, pmin, max, min)/2)
  expect_equivalent(y.summ[5], cyclicalweights(y.m, pmin, sum, min)/2)
  expect_equivalent(y.summ[6], cyclicalweights(y.m, pgeomean, sum, pgeomean)/2)

  # Check d_ statistics against the s_ statistics
  
  y.sim <- simulate.call(y)

  s_results <- t(sapply(y.sim, summary.call))
  d_results <- attr(y.sim,"stats")

  expect_equivalent(s_results,as.matrix(d_results))
})

test_that("valued triadic effects in directed networks", {
  
  # Directed
  y <- network.initialize(20, dir=TRUE)
  
  y <- simulate(y~sum,coef=0,reference=~DiscUnif(0,4),response="w",control=control.simulate(MCMC.burnin=1000),nsim=1)
  
  y.summ <- summary.call(y)
  
  y.m <- as.matrix(y, a="w")
  
  expect_equivalent(y.summ[1], transitiveweights(y.m, pmin, max, min))
  expect_equivalent(y.summ[2], transitiveweights(y.m, pmin, sum, min))
  expect_equivalent(y.summ[3], transitiveweights(y.m, pgeomean, sum, pgeomean))
  expect_equivalent(y.summ[4], cyclicalweights(y.m, pmin, max, min))
  expect_equivalent(y.summ[5], cyclicalweights(y.m, pmin, sum, min))
  expect_equivalent(y.summ[6], cyclicalweights(y.m, pgeomean, sum, pgeomean))
  
  # Check d_ statistics against the s_ statistics
  
  y.sim <- simulate.call(y)
  
  s_results <- t(sapply(y.sim, summary.call))
  d_results <- attr(y.sim,"stats")

  expect_equivalent(s_results,as.matrix(d_results))
})