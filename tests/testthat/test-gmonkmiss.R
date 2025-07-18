#  File tests/testthat/test-gmonkmiss.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

o <- options(ergm.eval.loglik=FALSE)

data(sampson)

run.test <- function() {
  #
  # Create random 25% missing
  #
  # msamplike <- rergm(network.size(samplike),prob=0.25, directed=FALSE)
  # samplike <- set.graph.attribute(samplike, "design", msamplike)
  # summary(samplike)

  set.seed(123)
  respondent <- rep(FALSE,network.size(samplike))
  respondent[sample(1:network.size(samplike), size=network.size(samplike)-2,replace=FALSE)] <- TRUE
  respondent

  summary(samplike)

  degreedist(samplike)

  set.seed(234)
  efit <- ergm(samplike~edges + gwesp(0.2, fixed=T), estimate="MPLE")
  summary(efit)

  set.seed(345)
  efit <- ergm(samplike~edges + gwesp(0.2, fixed=T), control=control.ergm(MCMLE.maxit=3))
  summary(efit)

  ## Test bounded degrees.
  set.seed(456)
  efit <- ergm(samplike~edges + gwesp(0.2, fixed=T), constraints=~bd(maxout=9),
               control=control.ergm(MCMLE.maxit=3))
  summary(efit)

  samplike <- set.vertex.attribute(samplike, "respondent", respondent)
  rm(respondent)
  summary(samplike)

  efit <- ergm(samplike~edges + gwesp(0.2, fixed=T), estimate="MPLE")
  summary(efit)

  set.seed(567)
  efit <- ergm(samplike~edges + gwesp(0.2, fixed=T), control=control.ergm(MCMLE.maxit=3))
  summary(efit)
}

test_that("directed network with missing data and dyadic dependence", {
  expect_error(run.test(), NA)
})

options(o)
