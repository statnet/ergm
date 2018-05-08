context("test-Sum.R")

data(florentine)

test_that("Sum() summary with one formula", {
  baseline <- summary(flomarriage~edges+absdiff("wealth"))
  test <- summary(flomarriage~Sum(cbind(1,.5)~edges+absdiff("wealth"),""))
  expect_equivalent(test, sum(baseline*c(1,.5)))
})


test_that("Sum() summary with one term and default weights", {
  baseline <- summary(flomarriage~absdiff("wealth"))
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), ~absdiff("wealth")),""))
  expect_equivalent(test, baseline*2)
})

test_that("Sum() summary with one term and differing weights", {
  baseline <- summary(flomarriage~absdiff("wealth"))
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), 0.5~absdiff("wealth")),""))
  expect_equivalent(test, baseline*1.5)
})


test_that("Sum() summary with default weights", {
  baseline <- summary(flomarriage~edges+absdiff("wealth"))
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")),""))
  expect_equivalent(test, baseline*2)
})

test_that("Sum() summary with differing weights", {
  baseline <- summary(flomarriage~edges+absdiff("wealth"))
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), 0.5~edges+absdiff("wealth")),""))
  expect_equivalent(test, baseline*1.5)
})

test_that("Sum() summary with heterogeneous lengths (error)", {
  expect_error(summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges),"")),"differ in length")
})

test_that("Sum() summary with matrix weights", {
  baseline <- summary(flomarriage~edges+absdiff("wealth"))
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), rbind(.5,0)~edges),""))
  expect_equivalent(test, baseline*c(1.5,1))
})

