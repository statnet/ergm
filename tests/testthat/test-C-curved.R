context("test-C-Curved.R")
data(faux.mesa.high)
flom <- ergm_model(~edges+gwesp()+gwdegree()+absdiffcat("Grade")+Offset(~nodefactor("Grade"),c(+1,-1), c(2,3))+gwesp()+NodematchFilter(~gwesp()+nodefactor("Grade"), "Grade"), faux.mesa.high)
(neta <- nparam(flom, canonical=TRUE))
(ntheta <- nparam(flom, canonical=FALSE))

test_that("C implementation of ergm.eta gives the same answer as R implementation.", {
  expect_equal(ergm.eta(1:ntheta, flom$etamap), ergm.eta.C(1:ntheta, flom$etamap))
})

test_that("C implementation of ergm.etagrad gives the same answer as R implementation.", {
  expect_equal(ergm.etagrad(1:ntheta, flom$etamap), ergm.etagrad.C(1:ntheta, flom$etamap))
})

test_that("C implementation of ergm.etagradmult gives the same answer as R implementation.", {
  v <- matrix(rnorm(2*neta), neta, 2)
  expect_equal(ergm.etagradmult(1:ntheta, v, flom$etamap), ergm.etagradmult.C(1:ntheta, v, flom$etamap))
})
