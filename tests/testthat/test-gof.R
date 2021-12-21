# Tests of gof()


test_that("gof() works on a bipartite ERGM (#424)", {
  mat <- matrix(c(1,1,0,1,1,0,1,0,0,1,0,1), 4, 3)
  net <-as.network(mat, bipartite=4)
  fit <- ergm(net ~ edges)
  expect_silent(
    r <- gof(fit)
  )
})
