#  File tests/testthat/test-fix.curved.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
test_that("fix.curved() works", {
  data(sampson)
  out<-fix.curved(samplike~edges+gwnsp(decay=.5,fixed=TRUE)+gwesp()+gwodegree()+edges,c(1:7))
  expect_equal(out,
               list(formula = samplike ~ edges + gwnsp(decay = 0.5, fixed = TRUE) + gwesp(fixed = TRUE, decay = 4L) + gwodegree(fixed = TRUE, decay = 6L) + edges,
                    theta = c(1,2,3,5,7)))
})
