expect_summary <- function(s, e, value, coefficients, tolerance=0.001) {
  expect_equal(s, value, tolerance=tolerance, ignore_attr=is.null(names(value)))
  expect_equal(coef(e)[1:length(coefficients)], coefficients, tolerance=tolerance, ignore_attr=is.null(names(coefficients)))
}
