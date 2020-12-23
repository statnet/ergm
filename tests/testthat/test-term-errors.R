test_that("ergm_Init_abort() can locate the term from which it had been called.", {
  data(faux.mesa.high)
  expect_error(summary(faux.mesa.high~nodefactor("Blah")), ".*term.*\\bnodefactor\\b.*in package.*\\bergm\\b.*")
})
