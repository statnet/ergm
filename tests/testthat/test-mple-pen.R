data(sampson)

test_that("Firth-penalized MPLE coefficients", {
  pen <- ergm(samplike ~ edges + nodefactor("group") + nodematch("group", diff = TRUE),
              control = control.ergm(MPLE.type = "penalized"))

  # library(logistf)
  # xy <- ergmMPLE(samplike ~ edges + nodefactor("group") + nodematch("group", diff = TRUE), output = "matrix")
  # firth <- logistf(formula = xy$response ~ xy$predictor - 1, weights = xy$weights)
  # firth |> coef() |> unname() |> deparse1() |> cat()
  expect_equal(coef(pen), c(-2.74802393589602, 0.0211052504894326,
                            0.985421050959423, 2.93460989147643,
                            4.14089796020648, 1.66917987328225),
               ignore_attr = TRUE, tolerance = 1e-6)
})
