context("test-bridge-target.stats.R")

library(ergm)
library(statnet.common)
logit <- function(p) log(p/(1-p))
expit <- function(e) exp(e)/(1+exp(e))
l <- function(nw, ets=NULL, theta=NULL){
  e <- NVL(ets, network.edgecount(nw))
  d <- network.dyadcount(nw)
  theta <- NVL(theta, logit(e/d))

  e*log(expit(theta)) + (d-e)*log(expit(-theta))
}

data(florentine)
y <- flomarriage

test_that("Log-likelihood with attainable target statistics",{
  set.seed(1)
  ts <- 3
  llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts)))
  llk <- l(y,ts)
  expect_equivalent(llk,llk.ergm,tolerance=0.01)
})

test_that("Log-likelihood with unattainable target statistics",{
  set.seed(1)
  ts <- 3.5
  llk.ergm <- as.vector(logLik(ergm(flomarriage~edges, target.stats=ts)))
  llk <- l(y,ts)
  expect_equivalent(llk,llk.ergm,tolerance=0.05)
})
