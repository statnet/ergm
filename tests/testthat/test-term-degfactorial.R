degm1lfactorial <- function(d){
  d <- d-1
  sum(lfactorial(d[d>=0]))
}

data(florentine)

test_that("degm1lfactorial summary", {
  expect_equal(summary(flomarriage~degm1lfactorial),
               degm1lfactorial(summary(flomarriage~sociality(nodes=TRUE))),
               ignore_attr = TRUE)
})

n <- 20
b <- 5
nw0 <- network.initialize(n, bipartite = b, directed = FALSE)
nw1 <- simulate(nw0 ~ edges, coef = 0)

test_that("b1degm1lfactorial summary", {
  expect_equal(summary(nw1~b1degm1lfactorial),
               degm1lfactorial(summary(nw1~b1sociality(nodes=TRUE))),
               ignore_attr = TRUE)
})

test_that("b2degm1lfactorial summary", {
  expect_equal(summary(nw1~b2degm1lfactorial),
               degm1lfactorial(summary(nw1~b2sociality(nodes=TRUE))),
               ignore_attr = TRUE)
})

data(sampson)

test_that("odegm1lfactorial summary", {
  expect_equal(summary(samplike~odegm1lfactorial),
               degm1lfactorial(summary(samplike~sender(nodes=TRUE))),
               ignore_attr = TRUE)
})

test_that("idegm1lfactorial summary", {
  expect_equal(summary(samplike~idegm1lfactorial),
               degm1lfactorial(summary(samplike~receiver(nodes=TRUE))),
               ignore_attr = TRUE)
})
