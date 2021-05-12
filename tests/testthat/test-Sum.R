
data(florentine)
baseline <- summary(flomarriage~edges+absdiff("wealth"))

test_that("Sum() summary with one formula", {
  test <- summary(flomarriage~Sum(cbind(1,.5)~edges+absdiff("wealth"),""))
  expect_equivalent(test, sum(baseline*c(1,.5)))
})


test_that("Sum() summary with one term and default weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), ~absdiff("wealth")),""))
  expect_equivalent(test, baseline[2]*2)
})

test_that("Sum() summary with one term and differing weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), 0.5~absdiff("wealth")),""))
  expect_equivalent(test, baseline[2]*1.5)
})


test_that("Sum() summary with default weights", {
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")),""))
  expect_equivalent(test, baseline*2)
})

test_that("Sum() summary with differing weights", {
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), 0.5~edges+absdiff("wealth")),""))
  expect_equivalent(test, baseline*1.5)
})

test_that("Sum() summary with heterogeneous lengths (error)", {
  expect_error(summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges),"")),"differ in length")
})

test_that("Sum() summary with matrix weights", {
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), rbind(.5,0)~edges),""))
  expect_equivalent(test, baseline*c(1.5,1))
})

test_that("Sum() summary with keyword weights", {
  test <- summary(flomarriage~Sum("sum"~edges+absdiff("wealth"),"")+Sum("mean"~edges+absdiff("wealth"),""))
  expect_equivalent(test, c(sum(baseline),mean(baseline)))
})


test_that("Prod() summary with default weights", {
  local_edition(3)
  test <- summary(flomarriage~Prod(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")),""))
  expect_equal(test, baseline^2, ignore_attr=TRUE)
})

test_that("Prod() summary with differing weights", {
  local_edition(3)
  test <- summary(flomarriage~Prod(c(~edges+absdiff("wealth"), 0.5~edges+absdiff("wealth")),""))
  expect_equal(test, baseline^1.5, ignore_attr=TRUE)
})

test_that("Prod() summary with keyword weights", {
  local_edition(3)
  test <- summary(flomarriage~Prod("prod"~edges+absdiff("wealth"),"")+Prod("geomean"~edges+absdiff("wealth"),""))
  expect_equal(test, c(prod(baseline),sqrt(prod(baseline))), ignore_attr=TRUE, tolerance=1e-5)
})
