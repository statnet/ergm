data(faux.mesa.high)
formula <- faux.mesa.high ~ nodematch("Sex")
mplesetup <- ergmMPLE(formula)

test_that("output = 'matrix' works for undirected networks", {
  local_edition(3)
  ord <- order(mplesetup$weights)
  m <- cbind(mplesetup$weights, mplesetup$response, mplesetup$predictor)[ord,]
  expect_equal(m, matrix(c(71, 132, 10284, 10423, 1, 1, 0, 0, 0, 1, 1, 0), 4,3), ignore_attr=TRUE)
})

test_that("output = 'fit' works for undirected networks", {
  local_edition(3)
  modelfit <- ergmMPLE(formula, output="fit")
  alt <- glm(mplesetup$response ~ mplesetup$predictor - 1,
             weights = mplesetup$weights, family="binomial")
  expect_equal(coef(modelfit), coef(alt), ignore_attr=TRUE, tolerance=0.001)
})
