context("test-c-ergm_model.R")

data(florentine)

test_that("concatenation of ergm_models works", {
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp(0,fixed=FALSE)+triangles+offset(gwdsp(0,fixed=FALSE))+absdiff("wealth"),flobusiness) -> m12
  ergm_model(~edges+offset(kstar(2))+absdiff("priorates")+gwesp(0,fixed=FALSE),flobusiness) -> m1
  ergm_model(~triangles+offset(gwdsp(0,fixed=FALSE))+absdiff("wealth"),flobusiness) -> m2
  expect_equal(c(m1,m2),m12)
})
