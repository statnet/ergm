
data(florentine)
baseline <- summary(flomarriage~edges+absdiff("wealth"))

test_that("Sum() summary with one formula, simple weights, procedural naming, and an offset", {
  local_edition(3)
  test <- summary(flomarriage~Sum(c(1,.5)~edges+offset(absdiff("wealth")), identity))
  expect_equal(test, setNames(baseline*c(1,.5), c("Sum~edges","offset(Sum~absdiff.wealth)")))
})

test_that("Sum() summary with one formula, simple weights, and procedural naming with AsIs", {
  local_edition(3)
  test <- summary(flomarriage~Sum(c(1,.5)~edges+absdiff("wealth"), I))
  expect_equal(test, setNames(baseline*c(1,.5), c("edges","absdiff.wealth")))
})

test_that("Sum() summary with one term and default weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), ~absdiff("wealth")),""))
  expect_equivalent(test, baseline[2]*2)
})

test_that("Sum() summary with one term and differing weights", {
  test <- summary(flomarriage~Sum(c(~absdiff("wealth"), 0.5~absdiff("wealth")),""))
  expect_equivalent(test, baseline[2]*1.5)
})

test_that("Sum() summary with default weights and procedural naming with AsIs", {
  local_edition(3)
  test <- summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges+absdiff("wealth")), function(x) I(paste(x[[1]], x[[2]]))))
  expect_equal(test, setNames(baseline*2, c("edges edges", "absdiff.wealth absdiff.wealth")))
})

test_that("Sum() summary with differing weights, offset (ignored), and forced-identical names", {
  expect_warning(test <- summary(flomarriage~Sum(c(~edges+offset(absdiff("wealth")), 0.5~edges+absdiff("wealth")),c("a","a"))), ".*does not propagate.*offset().*")
  expect_equivalent(test, setNames(baseline*1.5, c("a","a")))
})

test_that("Sum() summary with heterogeneous lengths (error)", {
  expect_error(summary(flomarriage~Sum(c(~edges+absdiff("wealth"), ~edges),"")),"differ in length")
})

test_that("Sum() summary with matrix weights and offset (ignored)", {
  local_edition(3)
  expect_warning(test <- summary(flomarriage~Sum(c(~edges+offset(absdiff("wealth")), rbind(.5,0)~edges),"")), ".*does not propagate.*offset().*")
  expect_equal(test, setNames(baseline*c(1.5,1), c("Sum~1", "Sum~2")))
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
