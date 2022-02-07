test_that("shrink_into_CH() works in both Rglpk and lpSolveAPI modes and produces identical results", {
  set.seed(0)
  p <- matrix(rnorm(200), 20)
  M <- matrix(rnorm(2000), 200)
  m <- colMeans(M)
  expect_equal(shrink_into_CH(p, M, m, solver="glpk"), shrink_into_CH(p, M, m, solver="lpsolve"))
})
