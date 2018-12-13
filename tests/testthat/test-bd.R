context("test-bd.R")

vary <- function(x, tol=sqrt(.Machine$double.eps)){
  apply(x,2,var)>tol
}

test_that("Bounded degree (bd()) maximum constraint for undirected networks", {
  data(florentine)
  constr <- simulate(flomarriage~edges, coef=1, monitor=~degree(0:network.size(flomarriage)), constraints=~bd(maxout=6), nsim=1000, output="stats")
  degr <- simulate(flomarriage~edges+degrange(7,Inf), coef=c(1,-Inf), monitor=~degree(0:network.size(flomarriage)), nsim=1000, output="stats")
  
  expect_true(approx.hotelling.diff.test(constr,degr[,-2])$p.value > 0.01)
  expect_true(all(vary(constr)==vary(degr[,-2])))
})

test_that("Bounded degree (bd()) constraints for directed networks", {
  data(sampson)
  constr <- simulate(samplike~edges, coef=-1, monitor=~odegree(0:network.size(samplike))+idegree(0:network.size(samplike)), constraints=~bd(minout=3,minin=2,maxout=6), nsim=1000, output="stats")
  degr <- simulate(samplike~edges+odegrange(0,3)+idegrange(0,2)+odegrange(7,Inf), coef=c(-1,-Inf, -Inf, -Inf), monitor=~odegree(0:network.size(samplike))+idegree(0:network.size(samplike)), nsim=1000, output="stats")
  
  expect_true(approx.hotelling.diff.test(constr,degr[,-(2:4)])$p.value > 0.01)
  expect_true(all(vary(constr)==vary(degr[,-(2:4)])))
})
