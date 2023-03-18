test_that("xTAx_qrssolve() works correctly in the presence of nullity", {
  x <- rnorm(2, sd=c(1,1e12))
  x <- c(x, sum(x))
  A <- matrix(c(1, 0, 1,
                0, 1e24, 1e24,
                1, 1e24, 1e24), 3, 3)

  expect_equal(xTAx_qrssolve(x,A),
               structure(drop(x%*%sginv(A)%*%x), rank = 2L, nullity = 1L))
})


test_that("xTAx_qrssolve() detects correctly if x is not in the span of A", {
  x <- rnorm(2, sd=c(1,1e12))
  x <- c(x, rnorm(1, sd=1e12))
  A <- matrix(c(1, 0, 1,
                0, 1e24, 1e24,
                1, 1e24, 1e24), 3, 3)

  expect_error(xTAx_qrssolve(x,A),
               "x is not in the span of A")
})
