context("test-c-ergm_model.R")

data(florentine)

test_that("concatenation of ergm_models works", {
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp(0,fixed=FALSE)+triangles+offset(gwdsp(0,fixed=FALSE))+absdiff("wealth"),flobusiness) -> m12
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp(0,fixed=FALSE),flobusiness) -> m1
  ergm_model(~triangles+offset(gwdsp(0,fixed=FALSE))+absdiff("wealth"),flobusiness) -> m2
  expect_equal(c(m1,m2),m12)
})

test_that("concatenation of ergm_models with redundant and nonredundant auxiliaries works", {
  ergm_model(~discord.inter.union.net(flomarriage), flomarriage) -> m1
  ergm_model(~discord.inter.union.net(flobusiness), flomarriage) -> m2
  ergm_model(~discord.inter.union.net(flomarriage), flomarriage) -> m3
  ergm_model(~discord.inter.union.net(flomarriage)+discord.inter.union.net(flobusiness)+discord.inter.union.net(flomarriage), flomarriage) -> m123
  expect_equal(c(m1,m2,m3),m123)
})
