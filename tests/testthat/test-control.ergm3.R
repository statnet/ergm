test_that("control.ergm3() defaults differ", {
  expect_equal(lapply(diff(control.ergm(), control.ergm3())[c("MCMLE.termination", "MCMLE.effectiveSize")], `[[`, "y"),
               list(MCMLE.termination = "Hummel", MCMLE.effectiveSize = NULL))
})

test_that("control.ergm3() is otherwise identical", {
  ## Two distinct tests are needed because the first one is a more
  ## natural use, but it affects the environment in which some
  ## parameters that are functions are defined.

  expect_equal(control.ergm(), control.ergm3(MCMLE.termination = "confidence", MCMLE.effectiveSize = 64), ignore_function_env = TRUE)

  expect_equal(control.ergm(), control.ergm3(MCMLE.termination = c("confidence", "Hummel", "Hotelling", "precision", "none"), MCMLE.effectiveSize = 64))
})
