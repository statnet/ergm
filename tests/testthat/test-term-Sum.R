#  File tests/testthat/test-term-Sum.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

data(florentine)
baseline <- summary(flomarriage~edges+absdiff("wealth"))
esps <- summary(flomarriage~gwesp(), gw.cutoff=4)

test_that("Sum() summary with one formula, simple weights, procedural naming, an offset, and a curved term", {
  test <- summary(flomarriage~Sum(c(1,.5, rep(2,4))~edges+offset(absdiff("wealth"))+gwesp(), identity), gw.cutoff=4)
  expect_equal(test, setNames(c(baseline,esps)*c(1,.5,rep(2,4)), c("Sum~edges","offset(Sum~absdiff.wealth)", "Sum~esp#1", "Sum~esp#2", "Sum~esp#3", "Sum~esp#4")))
})

test_that("Sum() summary with one formula, simple weights, fixed naming, an offset, and a curved term", {
  m <- ergm_model(flomarriage~Sum(c(1,.5, rep(2,4))~edges+offset(absdiff("wealth"))+gwesp(), "x"), gw.cutoff=4)
  expect_equal(param_names(m), c("Sum~x1", "offset(Sum~x2)", "Sum~x3", "Sum~x4"))
  expect_equal(summary(m, flomarriage), setNames(c(baseline,esps)*c(1,.5,rep(2,4)), replace(paste0("Sum~x",1:6), 2, "offset(Sum~x2)")))
})

test_that("Sum() summary with one formula, simple weights, fixed procedural naming with AsIs, an offset, and a curved term", {
  m <- ergm_model(flomarriage~Sum(c(1,.5, rep(2,4))~edges+offset(absdiff("wealth"))+gwesp(), I), gw.cutoff=4)
  expect_equal(param_names(m), c("edges", "offset(absdiff.wealth)", "gwesp", "gwesp.decay"))
  expect_equal(summary(m, flomarriage), setNames(c(baseline,esps)*c(1,.5,rep(2,4)), c("edges", "offset(absdiff.wealth)", "esp#1", "esp#2", "esp#3", "esp#4")))
})

test_that("Sum() summary with one formula, simple weights, and procedural naming with AsIs", {
  test <- summary(flomarriage~Sum(c(1,.5)~edges+absdiff("wealth"), I))
  expect_equal(test, setNames(baseline*c(1,.5), c("edges","absdiff.wealth")))
})

test_that("Sum() summary with one term and default weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), ~absdiff("wealth")),""))
  expect_equal(test, baseline[2]*2, ignore_attr=TRUE)
})

test_that("Sum() summary with one term and differing weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), 0.5~absdiff("wealth")),""))
  expect_equal(test, baseline[2]*1.5, ignore_attr=TRUE)
})

test_that("Sum() summary with default weights and procedural naming with AsIs", {
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")), function(x) I(paste(x[[1]], x[[2]]))))
  expect_equal(test, setNames(baseline*2, c("edges edges", "absdiff.wealth absdiff.wealth")))
})

test_that("Sum() summary with differing weights, offset (ignored), and forced-identical names", {
  expect_warning(test <- summary(flomarriage~Sum(c(~edges+offset(absdiff("wealth")), 0.5~edges+absdiff("wealth")),c("a","a"))), ".*does not propagate.*offset().*")
  expect_equal(test, setNames(baseline*1.5, c("a","a")), ignore_attr=TRUE)
})

test_that("Sum() summary with heterogeneous lengths (error)", {
  expect_error(summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges),"")),"differ in length")
})

test_that("Sum() summary with matrix weights and offset (ignored)", {
  expect_warning(test <- summary(flomarriage~Sum(c(~edges+offset(absdiff("wealth")), rbind(.5,0)~edges),"")), ".*does not propagate.*offset().*")
  expect_equal(test, setNames(baseline*c(1.5,1), c("Sum~1", "Sum~2")))
})

test_that("Sum() summary with keyword weights", {
  test <- summary(flomarriage~Sum("sum"~edges+absdiff("wealth"),"")+Sum("mean"~edges+absdiff("wealth"),""))
  expect_equal(test, c(sum(baseline),mean(baseline)), ignore_attr=TRUE)
})


test_that("Prod() summary with default weights", {
  test <- summary(flomarriage~Prod(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")),""))
  expect_equal(test, baseline^2, ignore_attr=TRUE)
})

test_that("Prod() summary with differing weights", {
  test <- summary(flomarriage~Prod(c(~edges+absdiff("wealth"), 0.5~edges+absdiff("wealth")),""))
  expect_equal(test, baseline^1.5, ignore_attr=TRUE)
})

test_that("Prod() summary with keyword weights", {
  test <- summary(flomarriage~Prod("prod"~edges+absdiff("wealth"),"")+Prod("geomean"~edges+absdiff("wealth"),""))
  expect_equal(test, c(prod(baseline),sqrt(prod(baseline))), ignore_attr=TRUE, tolerance=1e-5)
})
