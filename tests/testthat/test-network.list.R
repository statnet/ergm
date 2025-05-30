suppressMessages({
  net0 <- network.initialize(36)
  net0 %v% "x" <- rep(1:2, length = network.size(net0))
  k <- c(-2, -3, -3, -2)
  s1 <- simulate(net0 ~ mm("x", levels2=TRUE), nsim=10, coef = k)
  s2 <- simulate(net0 ~ mm("x", levels2=TRUE), nsim=10, coef = k * 0.01)
})

test_that("there is no error if attributes match", {
  expect_silent(
    c(s1, s1)
  )
  expect_silent(
    c(s1, s1)
  )
})

test_that("there is an error if attributes don't match", {
  expect_error(
    c(s1, s2)
  )
})

test_that("there is no error if attributes are to be ignored", {
  expect_silent(
    c(s1, s2, check_attr = FALSE)
  )
})
