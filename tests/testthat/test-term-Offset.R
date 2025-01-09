#  File tests/testthat/test-term-Offset.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
nw <- network(8, dir = FALSE)

test_that("Estimation with Offset() operator works", {
  offset_RE <- ".*offset.* decorator used on term .* with no free parameters is meaningless.*"

  expect_message(
    ergm_model(nw ~ edges + offset(Offset(~triangle, which = 1, coef = 0.1))),
    offset_RE)

  expect_failure(expect_message(
    off <- ergm(nw ~ edges + offset(triangle), offset.coef = 0.1),
    offset_RE), ".* did not throw the expected message.*")

  expect_failure(expect_message(
    Off <- ergm(nw ~ edges + Offset(~triangle, which = 1, coef = 0.1)),
    offset_RE), ".* did not throw the expected message.*")

  expect_lt(
  (coef(off)[1] - coef(Off)) /
  sqrt(vcov(off, sources = "estimation")[1,1] + vcov(Off, sources = "estimation")),
  qnorm(.9999))
})

test_that("Offset operator with curved terms works", {
  m <- ergm_model(nw ~ Offset(~edges + triangle + gwesp(0, fixed = T) + gwesp(fixed = F, cutoff = 3), which = 1, coef = 1))
  expect_equal(param_names(m), c("triangle", "gwesp.fixed.0", "gwesp", "gwesp.decay"))
  expect_equal(param_names(m, canonical=TRUE), c("edges", "triangle", "gwesp.fixed.0", "esp#1", "esp#2", "esp#3" ))
})
