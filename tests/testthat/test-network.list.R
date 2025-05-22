suppressMessages({
  data("faux.desert.high", package = "ergm")
  fit1 <- ergm(
    faux.desert.high ~ mm("sex"),
    control = control.ergm(
      force.main=TRUE,
      MCMLE.maxit = 1,
      seed = 1
    )
  )
  fit2 <- ergm(
    faux.desert.high ~ mm("sex"),
    control = control.ergm(
      force.main=TRUE,
      MCMLE.maxit = 1,
      seed = 2
    )
  )
  s1 <- simulate(fit1, 10)
  s2 <- simulate(fit2, 10)
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
