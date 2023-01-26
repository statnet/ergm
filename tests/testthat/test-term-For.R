#  File tests/testthat/test-term-For.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
test_that("For() operator with list input", {
  data(florentine)
  expect_equal(summary(flomarriage ~ For(~absdiff("wealth", pow=.), . = 1:3)),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3)))
})

test_that("For() operator with formula input", {
  data(sampson)
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = ~sort(unique(.%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
})

test_that("For() operator with function input", {
  data(sampson)
  expect_equal(summary(samplike ~ For(~nodematch("group", levels=., diff=TRUE), . = function(nw) sort(unique(nw%v%"group")))),
               summary(samplike ~ nodematch("group", diff=TRUE)))
})

test_that("For() operator with multiple list inputs", {
  data(florentine)
  expect_equal(summary(flomarriage ~ For(~absdiff(a, pow=.), a = c("wealth", "priorates"), . = 1:3)),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("wealth", pow=2) + absdiff("wealth", pow=3) +
                         absdiff("priorates", pow=1) + absdiff("priorates", pow=2) + absdiff("priorates", pow=3)))

  expect_equal(summary(flomarriage ~ For(. = 1:3, a = c("wealth", "priorates"), ~absdiff(a, pow=.))),
               summary(flomarriage ~ absdiff("wealth", pow=1) + absdiff("priorates", pow=1) +
                         absdiff("wealth", pow=2) + absdiff("priorates", pow=2) +
                         absdiff("wealth", pow=3) + absdiff("priorates", pow=3)))
})
