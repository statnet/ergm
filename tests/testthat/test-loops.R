o <- options(ergm.warn_loops = FALSE)

test_that("summary for degree and sociality with loops", {
  nw0 <- network.initialize(10, loops = TRUE, directed = FALSE)
  set.seed(0)
  nw1 <- simulate(nw0 ~ edges, coef = -1)
  degs <- rowSums(as.matrix(nw1))

  expect_equal(summary(nw1 ~ degree(0:10)), tabulate(degs + 1, 11),
               ignore_attr = TRUE)
  expect_equal(summary(nw1 ~ sociality(nodes = TRUE)), degs,
               ignore_attr = TRUE)
})

options(o)
