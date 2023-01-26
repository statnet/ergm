#  File tests/testthat/test-C-curved.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

ergm.eta.R <- function(theta, etamap) {
  eta <- numeric(etamap$etalength)
  ec <- etamap$canonical
  eta[ec[ec>0]] <- theta[ec>0]
  if(length(etamap$curved)>0) {
    for(cm in etamap$curved) {
      eta[cm$to] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
    }
  }
  eta
}

ergm.etagrad.R <- function(theta, etamap) {
  etagrad <- matrix(0, length(theta), etamap$etalength)
  ec <- etamap$canonical
# Set gradient for canonical parameters to the identity matrix
  etagrad[ec>0, ec[ec>0]] <- diag(sum(ec>0))
  if(length(etamap$curved)>0) {
    for(cm in etamap$curved) {
      etagrad[cm$from,cm$to] <- cm$gradient(theta[cm$from], length(cm$to), cm$cov)
    }
  }
  etagrad
}

ergm.etagradmult.R <- function(theta, v, etamap) {
  v <- as.matrix(v)
  ans <- matrix(0, length(theta), dim(v)[2])
  if(dim(v)[1] != etamap$etalength)
    stop("Non-conforming matrix multiply: grad(eta) %*% v.\n",
         "grad(eta) has ", etamap$etalength, " columns ",
         "and v has ", dim(v)[1], " rows.")
  ec <- etamap$canonical
# Set gradient for canonical parameters to the identity matrix
  ans[ec>0,] <- v[ec[ec>0],]
  if(length(etamap$curved)>0) {
    for(cm in etamap$curved) {
      ans[cm$from,] <- cm$gradient(theta[cm$from], length(cm$to), cm$cov)%*%v[cm$to,]
    }
  }
  ans
}

data(faux.mesa.high)
flom <- ergm_model(~edges+gwesp()+gwdegree()+absdiffcat("Grade")+Offset(~nodefactor("Grade"),c(+1,-1), c(2,3))+gwesp()+NodematchFilter(~gwesp()+nodefactor("Grade"), "Grade"), faux.mesa.high)
(neta <- nparam(flom, canonical=TRUE))
(ntheta <- nparam(flom, canonical=FALSE))

test_that("C implementation of ergm.eta gives the same answer as R implementation.", {
  expect_equal(ergm.eta(1:ntheta, flom$etamap), ergm.eta.R(1:ntheta, flom$etamap))
})

test_that("C implementation of ergm.etagrad gives the same answer as R implementation.", {
  expect_equal(ergm.etagrad(1:ntheta, flom$etamap), ergm.etagrad.R(1:ntheta, flom$etamap))
})

test_that("C implementation of ergm.etagradmult gives the same answer as R implementation.", {
  v <- matrix(rnorm(2*neta), neta, 2)
  expect_equal(ergm.etagradmult(1:ntheta, v, flom$etamap), ergm.etagradmult.R(1:ntheta, v, flom$etamap))
})
