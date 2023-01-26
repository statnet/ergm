#  File tests/testthat/test-bd.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

set.seed(0)

vary <- function(x, tol=sqrt(.Machine$double.eps)){
  apply(x,2,function(x)length(unique(x)))
}

compvary <- function(x, y, diff.tol = 2, vary.tol=sqrt(.Machine$double.eps)){
  x <- vary(x, vary.tol)
  y <- vary(y, vary.tol)

  (x-y)^2/(x+y)*2 <= diff.tol
}

test_that("Bounded degree (bd()) maximum constraint for undirected networks", {
  data(florentine)
  constr <- simulate(flomarriage~edges, coef=0.5, monitor=~degree(0:network.size(flomarriage)), constraints=~bd(maxout=6), nsim=200, output="stats")
  degr <- simulate(flomarriage~edges+degrange(7,Inf), coef=c(0.5,-Inf), monitor=~degree(0:network.size(flomarriage)), nsim=200, output="stats")
  
  expect_warning(expect_gt(approx.hotelling.diff.test(constr,degr[,-2])$p.value, 0.01), ".*do not vary but equal mu0.*")
  expect_true(all(compvary(constr,degr[,-2])))
})

test_that("Bounded degree (bd()) constraints for directed networks", {
  data(sampson)
  constr <- simulate(samplike~edges, coef=-1, monitor=~odegree(0:network.size(samplike))+idegree(0:network.size(samplike)), constraints=~bd(minout=3,minin=2,maxout=6), nsim=200, output="stats")
  degr <- simulate(samplike~edges+odegrange(0,3)+idegrange(0,2)+odegrange(7,Inf), coef=c(-1,-Inf, -Inf, -Inf), monitor=~odegree(0:network.size(samplike))+idegree(0:network.size(samplike)), nsim=200, output="stats")
  
  expect_warning(expect_gt(approx.hotelling.diff.test(constr,degr[,-(2:4)])$p.value, 0.01), ".*do not vary but equal mu0.*")
  expect_true(all(compvary(constr,degr[,-(2:4)])))
})
