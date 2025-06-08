#  File tests/testthat/test-c-ergm_model.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

data(florentine)

test_that("concatenation of ergm_models works", {
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp()+triangles+offset(gwdsp())+absdiff("wealth"),flobusiness) -> m12
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp(),flobusiness) -> m1
  ergm_model(~triangles+offset(gwdsp())+absdiff("wealth"),flobusiness) -> m2
  cm12 <- c(m1,m2)
  cm12$uid <- m12$uid <- NULL
  expect_equal(cm12, m12, ignore_formula_env=TRUE)
})

test_that("concatenation of ergm_models with redundant and nonredundant auxiliaries works", {
  ergm_model(~discord.inter.union.net(flomarriage), flomarriage) -> m1
  ergm_model(~discord.inter.union.net(flobusiness), flomarriage) -> m2
  ergm_model(~discord.inter.union.net(flomarriage), flomarriage) -> m3
  ergm_model(~discord.inter.union.net(flomarriage)+discord.inter.union.net(flobusiness)+discord.inter.union.net(flomarriage), flomarriage) -> m123
  cm123 <- c(m1,m2,m3)
  cm123$uid <- m123$uid <- NULL
  expect_equal(cm123, m123, ignore_formula_env=TRUE)
})
