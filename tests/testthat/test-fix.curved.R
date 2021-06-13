test_that("fix.curved() works", {
  local_edition(3)
  data(sampson)
  out<-fix.curved(samplike~edges+gwnsp(decay=.5,fixed=TRUE)+gwesp()+gwodegree()+edges,c(1:7))
  expect_equal(out,
               list(formula = samplike ~ edges + gwnsp(decay = 0.5, fixed = TRUE) + gwesp(fixed = TRUE, decay = 4L) + gwodegree(fixed = TRUE, decay = 6L) + edges,
                    theta = c(1,2,3,5,7)))
})
